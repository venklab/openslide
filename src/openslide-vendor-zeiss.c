/*
 *  OpenSlide, a library for reading whole slide image files
 *
 *  Copyright (c) 2007-2013 Carnegie Mellon University
 *  Copyright (c) 2011 Google, Inc.
 *  All rights reserved.
 *
 *  OpenSlide is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, version 2.1.
 *
 *  OpenSlide is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with OpenSlide. If not, see
 *  <http://www.gnu.org/licenses/>.
 *
 * Add Zeiss czi support. Wei Chen <chenw1@uthscsa.edu>
 *
 */
#include <config.h>

#include "openslide-private.h"
#include "openslide-decode-jp2k.h"
#include "openslide-decode-tiff.h"
#include "openslide-decode-tifflike.h"
#include "openslide-decode-xml.h"
#include "openslide-decode-jxr.h"

#include <glib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <tiffio.h>
#include <libxml/tree.h>
#include <libxml/xpath.h>

/* from linux kernel math.h */
#define DIV_ROUND_CLOSEST(x, divisor)(             \
{                                                  \
        typeof(x) __x = x;                         \
        typeof(divisor) __d = divisor;             \
        (((typeof(x))-1) > 0 ||                    \
         ((typeof(divisor))-1) > 0 ||              \
         (((__x) > 0) == ((__d) > 0))) ?           \
                (((__x) + ((__d) / 2)) / (__d)) :  \
                (((__x) - ((__d) / 2)) / (__d));   \
}                                                  \
)

#define CZI_SEG_ID_LEN   16

//static const char CZI_ZISRAW[] = "ZISRAW";

/* zeiss uses little-endian */
struct zisraw_seg_hdr {
  char sid[CZI_SEG_ID_LEN];
  int64_t allocated_size;
  int64_t used_size;
};

struct __attribute__ ((__packed__)) zisraw_data_file_hdr {
  int32_t major;
  int32_t minor;
  int32_t _reserved1;
  int32_t _reserved2;
  char primary_file_guid[16];
  char file_guid[16];
  int32_t file_part;   // this causes off-align
  int64_t subblk_dir_pos;
  int64_t meta_pos;
  int32_t update_pending;
  int64_t attach_dir_pos;
};

struct zisraw_seg_subblk_hdr {
  struct zisraw_seg_hdr seg_hdr;
  int32_t entry_count;
  char _reserved[124];
};

struct zisraw_seg_meta_hdr {
  struct zisraw_seg_hdr seg_hdr;
  int32_t xml_size;
  int32_t _attach_size;
  char _reserved[248];
};

struct zisraw_data_subblk_hdr {
  int32_t meta_size;
  int32_t attach_size;
  int64_t data_size;
};

struct __attribute__ ((__packed__)) zisraw_dir_entry_dv {
  char schema[2];
  int32_t pixel_type;
  int64_t file_pos;
  int32_t _file_part;
  int32_t compression;
  int8_t pyramid_type;
  char _reserved1;
  char _reserved2[4];
  int32_t ndimensions;
  // followed by variable length array of zisraw_dim_entry_dv
};

struct zisraw_dim_entry_dv {
  char dimension[4];
  int32_t start;
  int32_t size;
  float start_coordinate;
  int32_t stored_size;
};

struct czi_subblk {
  int64_t file_pos;
  int64_t downsample_i;
  int32_t pixel_type;
  int32_t compression;
  int32_t x1, x2, y1, y2;
  uint32_t w, h, tw, th;
  int32_t id;
  int32_t dir_entry_len;
  int8_t channel;
  int8_t scene;
};

enum z_pyramid_type {
  PYR_NONE = 0,
  PYR_SINGLE,
  PYR_MULTIPLE,
};

enum z_compression {
  COMP_NONE = 0,
  COMP_JPEG,
  COMP_LZW,
  COMP_JXR = 4,
  COMP_OTHER,
};

enum pixel_type {
  PT_GRAY8 = 0,
  PT_GRAY16,
  PT_GRAY32FLOAT,
  PT_BGR24,
  PT_BGR96FLOAT = 8,
  PT_BGRA32,
  PT_GRAY64COMPLEX,
  PT_BGR192COMPLEX,
  PT_GRAY32,
  PT_GRAY64,
};

struct zeiss_ops_data {
  GMutex mutex;
  char *filename;
  int64_t subblk_dir_pos;
  int64_t meta_pos;
  int64_t attach_dir_pos;
  int32_t nsubblk;  // total number of subblocks
  int32_t w;
  int32_t h;
  int32_t scene;
  int32_t offset_x;
  int32_t offset_y;
  GPtrArray *subblks;
  GHashTable *grids;

  GHashTable *count_tile_width;
  GHashTable *count_tile_height;
  void *unused;
};

struct freq_count {
  int64_t value;
  int64_t count;
};

/*
static struct zeiss_ops_data zpd = {
  .filename = NULL,
  .subblk_dir_pos = 0,
  .meta_pos = 0,
  .attach_dir_pos = 0,
  .nsubblk = 0,
  .offset_x = G_MAXINT32,
  .offset_y = G_MAXINT32,
};
*/

struct level {
  struct _openslide_level base;
  int64_t downsample_i;
  int32_t compression;
};

static void destroy(openslide_t *osr);
static bool paint_region(openslide_t *osr, cairo_t *cr, int64_t x, int64_t y,
                         struct _openslide_level *level, int32_t w, int32_t h,
                         GError **err);

static const struct _openslide_ops zeiss_ops = {
  .paint_region = paint_region,
  .destroy = destroy,
};

static void destroy_level(struct level *l) {
  g_slice_free(struct level, l);
}

static void destroy_subblk(struct czi_subblk *p)
{
  g_slice_free(struct czi_subblk, p);
}

static void destroy(openslide_t *osr) {
  for (int32_t i = 0; i < osr->level_count; i++) {
    destroy_level((struct level *) osr->levels[i]);
  }
  g_free(osr->levels);

  struct zeiss_ops_data *data = osr->data;
  g_free(data->filename);
  g_hash_table_destroy(data->count_tile_width);
  g_hash_table_destroy(data->count_tile_height);
  g_hash_table_destroy(data->grids);

  g_ptr_array_free(data->subblks, TRUE);
  g_slice_free(struct zeiss_ops_data, data);
}

/* locate a grid based on level and let grid_paint_region do the job */
static bool paint_region(openslide_t *osr, cairo_t *cr,
                          int64_t x, int64_t y,
                          struct _openslide_level *level,
                          int32_t w, int32_t h,
                          GError **err)
{
  //int32_t w G_GNUC_UNUSED, int32_t h G_GNUC_UNUSED,
  struct zeiss_ops_data *data = osr->data;
  struct level *l = (struct level *) level;

  int64_t ds = (int64_t) level->downsample;
  printf("debug paint_region: ds = %ld, x, y (%ld, %ld), w, h %d/%d\n",
         ds, x, y, w, h);
  struct _openslide_grid *grid = g_hash_table_lookup(data->grids, &ds);
  void *unused_args = NULL;

  // need convert level 0 x,y to x,y on current level
  if (!_openslide_grid_paint_region(grid, cr, &unused_args,
                                    x / l->base.downsample,
                                    y / l->base.downsample,
                                    level, w, h, err)) {
      return false;
  }

  return true;
}


static bool zeiss_detect(const char *filename,
                         struct _openslide_tifflike *tl, GError **err) {
  // reject TIFFs
  if (tl) {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED,
                "Is a TIFF file");
    return false;
  }

  g_autoptr(_openslide_file) f = _openslide_fopen(filename, err);
  if (!f)
    return false;

  // string ZISRAWFILE occurs once per file, at positon 0
  char sid[CZI_SEG_ID_LEN];
  if (_openslide_fread(f, sid, CZI_SEG_ID_LEN) != CZI_SEG_ID_LEN) {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED,
                "Can't read CZI file header");
    return false;
  }

  if (!g_str_has_prefix(sid, "ZISRAWFILE")) {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED,
                "Not a Zeiss CZI slide");
    return false;
  }

  //TODO check pyramid

  return true;
}

static void read_dim_entry(struct zisraw_dim_entry_dv *p,
                           struct czi_subblk *sb)
{
  int start = GINT32_FROM_LE(p->start);
  int size = GINT32_FROM_LE(p->size);
  int stored_size = GINT32_FROM_LE(p->stored_size);

  switch (p->dimension[0]) {
  case 'X':
    sb->x1 = start;
    sb->w = size;
    sb->tw = stored_size;
    sb->x2 = start + size - 1;
    sb->downsample_i = DIV_ROUND_CLOSEST(size, stored_size);
    break;

  case 'Y':
    sb->y1 = start;
    sb->h = size;
    sb->th = stored_size;
    sb->y2 = start + size - 1;
    break;

  case 'C':
    sb->channel = start;
    break;

  case 'S':
    sb->scene = start;
    break;
  }
}

static void print_tile(struct czi_subblk *tile)
{
  printf("debug tile: s%dch%d, ds%ld, x1,y1(%d, %d), w/tw %d/%d, h/th %d/%d, pos %ld\n",
      tile->scene,
      tile->channel,
      tile->downsample_i,
      tile->x1,
      tile->y1,
      tile->w,
      tile->tw,
      tile->h,
      tile->th,
      tile->file_pos
      );

}

static int read_dir_entry(GPtrArray *subblks, char *p)
{
  struct czi_subblk *sb = g_slice_new(struct czi_subblk);
  char *b = p;

  struct zisraw_dir_entry_dv *dv = (struct zisraw_dir_entry_dv *) b;
  sb->pixel_type = GINT32_FROM_LE(dv->pixel_type);
  sb->compression = GINT32_FROM_LE(dv->compression);
  sb->file_pos = GINT64_FROM_LE(dv->file_pos);

  int nread = sizeof(struct zisraw_dir_entry_dv);
  int32_t ndim = GINT32_FROM_LE(dv->ndimensions);
  b += nread;  // the first entry of dimensions
  for (int i = 0; i < ndim; i++) {
    read_dim_entry((struct zisraw_dim_entry_dv *)b , sb);
    b += 20;
  }

  //print_blk(&blk);
  nread += ndim * 20;
  sb->dir_entry_len = nread;
  g_ptr_array_add(subblks, sb);

  return nread;
}


static bool read_subblk_dir(openslide_t *osr, GError **err)
{
  struct zeiss_ops_data *data = osr->data;
  data->subblks = g_ptr_array_new_full(64, (GDestroyNotify) destroy_subblk);

  g_autoptr(_openslide_file) f = _openslide_fopen(data->filename, err);
  if (!f)
    return false;

  char buf[512];
  if (!_openslide_fseek(f, data->subblk_dir_pos, SEEK_SET, err)) {
    g_prefix_error(err, "Couldn't seek to SubBlock directory");
    return false;
  }

  size_t len;
  len  = sizeof(struct zisraw_seg_subblk_hdr);
  if (_openslide_fread(f, buf, len) != len) {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED,
                "Cannot read FileHeader");
    return false;
  }

  struct zisraw_seg_subblk_hdr *hdr = (struct zisraw_seg_subblk_hdr *) buf;
  if (!g_str_has_prefix(hdr->seg_hdr.sid, "ZISRAWDIRECTORY")) {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED,
                "Not SubBlockDirectory");
    return false;
  }

  data->nsubblk = GINT32_FROM_LE(hdr->entry_count);
  printf("debug total %d subblocks\n", data->nsubblk);

  len = (size_t) GINT32_FROM_LE(hdr->seg_hdr.allocated_size);
  len -= sizeof(struct zisraw_seg_subblk_hdr);
  g_autofree char *buf_dir = g_slice_alloc(len);
  if (_openslide_fread(f, buf_dir, len) != len) {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED,
                "Cannot read SubBlockDirectory");
    return false;
  }

  char *p = buf_dir;
  int dir_entry_len;
  size_t total = 0;
  for (int i = 0; i < data->nsubblk; i++) {
    if (total > len) {
      g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED,
                  "Read beyond last byte of directory entires");
      return false;
    }

    dir_entry_len = read_dir_entry(data->subblks, p);
    p += dir_entry_len;
    total += dir_entry_len;
  }

  return true;
}

/* the topleft-most tile has none-zero (x, y), use its x,y as offset to adjust
 * x,y of other tiles.
 */
static bool adjust_coordinate_origin(openslide_t *osr,
                                     GError **err G_GNUC_UNUSED)
{
  struct zeiss_ops_data *data = osr->data;
  GPtrArray *subblks = data->subblks;
  struct czi_subblk *b;

  for (guint i = 0; i < subblks->len; i++) {
    b = subblks->pdata[i];
    if (b->x1 < data->offset_x)
      data->offset_x = b->x1;

    if (b->y1 < data->offset_y)
      data->offset_y = b->y1;
  }

  printf("debug: set offset x,y = %d, %d\n", data->offset_x, data->offset_y);
  for (guint i = 0; i < subblks->len; i++) {
    b = subblks->pdata[i];
    b->x1 -= data->offset_x;
    b->y1 -= data->offset_y;
  }

  return true;
}

static bool read_data_from_subblk(openslide_t *osr, struct czi_subblk *sb,
                                  void *dest, GError **err)
{
  // only work with BGR24 for now
  if (sb->pixel_type != PT_BGR24) {
    printf("debug read_data_from_subblk: wrong pixel_type\n");
    return false;
  }

  //printf("debug: subblock file_pos = %ld\n", sb->file_pos);
  struct zeiss_ops_data *data = osr->data;
  g_autoptr(_openslide_file) f = _openslide_fopen(data->filename, err);
  if (!f)
    return false;

  int64_t pos = sb->file_pos + 32;
  if (!_openslide_fseek(f, pos, SEEK_SET, err)) {
    g_prefix_error(err, "Couldn't seek to SubBlock");
    return false;
  }

  char buf[512];
  size_t len;
  len  = sizeof(struct zisraw_data_subblk_hdr);
  if (_openslide_fread(f, buf, len) != len)
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED,
                "Cannot read SubBlock header");

  struct zisraw_data_subblk_hdr *hdr = (struct zisraw_data_subblk_hdr *)buf;
  int n = MAX(256 - 16 - sb->dir_entry_len, 0);
  int offset_meta = MAX(256, n);
  pos += offset_meta + GINT32_FROM_LE(hdr->meta_size);
  //printf("debug: datalen = %ld, datastart = %ld\n", hdr->data_size,pos);
  if (!_openslide_fseek(f, pos, SEEK_SET, err)) {
    g_prefix_error(err, "Couldn't seek to pixel data");
    return false;
  }

  len = (size_t) GINT64_FROM_LE(hdr->data_size);
  //printf("debug: subblock meta len = %d, data len = %d\n",
  //    hdr->meta_size, hdr->data_size);
  g_autofree char *buf2 = (char *) g_slice_alloc(len);
  if (_openslide_fread(f, buf2, len) != len) {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED,
                "Cannot read pixel data");
    return false;
  }

  switch (sb->compression) {
  case COMP_JXR:
    _openslide_jxr_decode_buf(buf2, len, (struct decoded_jxr *) dest, NULL);
    break;

  case COMP_JPEG:
    break;

  case COMP_LZW:
    break;

  default:
    break;
  }

  return true;
}

static bool read_tile(openslide_t *osr, cairo_t *cr,
                      struct _openslide_level *level G_GNUC_UNUSED,
                      int64_t tid G_GNUC_UNUSED, void *data,
                      void *arg G_GNUC_UNUSED, GError **err G_GNUC_UNUSED)
{
  struct czi_subblk *sb = data;

  /*
  fprintf(stderr, "debug read_tile: pos = %ld, ds %ld, x1,y1(%d,%d), tw,th = %d, %d\n",
      sb->file_pos, sb->downsample_i, sb->x1, sb->y1, sb->tw, sb->th);
      */

  g_autoptr(_openslide_cache_entry) cache_entry = NULL;
  // file_pos is unique
  uint8_t *img = (uint8_t *) _openslide_cache_get(osr->cache, 0, sb->file_pos, 0,
                                                  &cache_entry);
  struct decoded_jxr dest;
  if (!img) {
    fprintf(stderr, "cache missing\n");
    read_data_from_subblk(osr, sb, &dest, NULL);
    g_assert(sb->tw == dest.w);
    g_assert(sb->th == dest.h);

    // _openslide_cache_entry_unref will free data
    _openslide_cache_put(osr->cache, 0, sb->file_pos, 0, dest.data, dest.size,
                         &cache_entry);
    img = dest.data;
  }

  int stride = cairo_format_stride_for_width(CAIRO_FORMAT_RGB24, sb->tw);

  g_autoptr(cairo_surface_t) surface =
    cairo_image_surface_create_for_data((unsigned char *) img,
                                        CAIRO_FORMAT_RGB24,
                                        sb->tw, sb->th, stride);
  cairo_set_source_surface(cr, surface, 0, 0);
  //cairo_surface_destroy(surface);
  cairo_paint(cr);

  return true;
}

static void finish_adding_tiles(void *key G_GNUC_UNUSED, void *value,
                                void *user_data G_GNUC_UNUSED)
{
  struct _openslide_grid *grid = (struct _openslide_grid *) value;
  _openslide_grid_range_finish_adding_tiles(grid);
}

static void destroy_int64_key(void *p)
{
  g_slice_free(int64_t, p);
}

static void destroy_freq_count(void *p)
{
  g_slice_free(struct freq_count, p);
}

static void destroy_wh_count_hashtable(void *p)
{
  g_hash_table_destroy((GHashTable *) p);
}

/* count occurrence of tile width and height for each level */
static void count_tile_width_height(openslide_t *osr, int64_t downsample,
                                    int64_t w, int64_t h)
{
  struct zeiss_ops_data *data = osr->data;
  struct freq_count *frq_w, *frq_h;
  GHashTable *level_w, *level_h;
  int64_t *k;

  // one hashtable for each level
  level_w = g_hash_table_lookup(data->count_tile_width, &downsample);
  if (!level_w) {
    level_w = g_hash_table_new_full(g_int64_hash, g_int64_equal,
                                    (GDestroyNotify) destroy_int64_key,
                                    (GDestroyNotify) destroy_freq_count);
    k = g_new(int64_t, 1);
    *k = downsample;
    g_hash_table_insert(data->count_tile_width, k, level_w);
    //fprintf(stderr, "created count_tile_width ht for ds %ld\n", downsample);
  }

  level_h = g_hash_table_lookup(data->count_tile_height, &downsample);
  if (!level_h) {
    level_h = g_hash_table_new_full(g_int64_hash, g_int64_equal,
                                    (GDestroyNotify) destroy_int64_key,
                                    (GDestroyNotify) destroy_freq_count);
    k = g_new(int64_t, 1);
    *k = downsample;
    g_hash_table_insert(data->count_tile_height, k, level_h);
    //fprintf(stderr, "created count_tile_height ht for ds %ld\n", downsample);
  }

  // one struct freq_count for each tile width in this level
  frq_w = g_hash_table_lookup(level_w, &w);
  if (!frq_w) {
    frq_w = g_slice_new0(struct freq_count);
    k = g_new(int64_t, 1);
    *k = w;
    g_hash_table_insert(level_w, k, frq_w);
    frq_w->value = w;
  }

  frq_h = g_hash_table_lookup(level_h, &h);
  if (!frq_h) {
    frq_h = g_slice_new0(struct freq_count);
    k = g_new(int64_t, 1);
    *k = h;
    g_hash_table_insert(level_h, k, frq_h);
    frq_h->value = h;
  }

  frq_w->count++;
  frq_h->count++;
}


static bool init_range_grids(openslide_t *osr, GError **err G_GNUC_UNUSED)
{
  struct zeiss_ops_data *data = osr->data;
  GPtrArray *subblks = data->subblks;
  struct czi_subblk *b;
  struct _openslide_grid *grid;
  struct level *l;
  int64_t *k;

  data->grids = g_hash_table_new_full(g_int64_hash, g_int64_equal,
                                      (GDestroyNotify) destroy_int64_key,
                                      (GDestroyNotify) _openslide_grid_destroy);
  for (int i = 0; i < osr->level_count; i++) {
    l = (struct level *) osr->levels[i];
    fprintf(stderr, "debug: level %d, ds %ld, common tw,th = %ld,%ld\n",
        i, l->downsample_i, l->base.tile_w, l->base.tile_h);
    grid = _openslide_grid_create_range(osr, l->base.tile_w, l->base.tile_h,
                                        read_tile,
                                        NULL);
    k = g_new(int64_t, 1);
    *k = l->downsample_i;
    g_hash_table_insert(data->grids, k, grid);
  }

  for (guint i = 0; i < subblks->len; i++) {
    b = subblks->pdata[i];
    grid = g_hash_table_lookup(data->grids, &b->downsample_i);
    _openslide_grid_range_add_tile(grid,
       (double) b->x1 / (double) b->downsample_i ,
       (double) b->y1 / (double) b->downsample_i ,
       (double) b->tw, (double) b->th, b);
  }

  fprintf(stderr, "deubg: finished init grids\n");
  g_hash_table_foreach(data->grids, (GHFunc)finish_adding_tiles, NULL);
  return true;
}

static gint cmp_int64(gpointer a, gpointer b) {
  int64_t *x = (int64_t *)a;
  int64_t *y = (int64_t *)b;
  if (*x == *y)
    return 0;

  return (*x < *y) ?  -1 : 1;
}

static void iter_common(void *k G_GNUC_UNUSED, void *value, void *user_data)
{
  struct freq_count *frq = (struct freq_count *) value;
  struct freq_count *result = (struct freq_count *) user_data;

  if (frq->count > result->count) {
    result->count = frq->count;
    result->value = frq->value;
  }
}

static inline int64_t find_most_common(GHashTable *ht, int64_t downsample)
{
  GHashTable *level;
  struct freq_count cnt = {.count = 0};

  //fprintf(stderr, "find_most_common: ds %ld level ht has %d keys\n",
  //    downsample, g_hash_table_size(ht));
  level = g_hash_table_lookup(ht, &downsample);
  if (!level) {
    fprintf(stderr, "find_most_common: ds %ld level is NULL\n", downsample);
    return 1;
  }

  g_hash_table_foreach(level, (GHFunc)iter_common, &cnt);
  //printf("debug: tw th %ld has %ld tiles\n", cnt.value, cnt.count);
  return cnt.value;
}

static int64_t find_most_common_width(openslide_t *osr, int64_t downsample)
{
  struct zeiss_ops_data *data = osr->data;
  return find_most_common(data->count_tile_width, downsample);
}

static int64_t find_most_common_height(openslide_t *osr, int64_t downsample)
{
  struct zeiss_ops_data *data = osr->data;
  return find_most_common(data->count_tile_height, downsample);
}

static bool init_levels(openslide_t *osr, GError **err G_GNUC_UNUSED)
{
  struct zeiss_ops_data *data = osr->data;
  GPtrArray *subblks = data->subblks;
  struct czi_subblk *b;
  struct level *l;

  for (guint i = 0; i < subblks->len; i++) {
    b = subblks->pdata[i];
    count_tile_width_height(osr, b->downsample_i, b->tw, b->th);
  }

  GList *downsamples = g_hash_table_get_keys(data->count_tile_width);
  GList *p = g_list_sort(downsamples, (GCompareFunc) cmp_int64);
  downsamples = p;

  g_autoptr(GPtrArray) levels = g_ptr_array_new();
  int64_t downsample_i;
  while (p) {
    //printf("debug: downsample_i = %ld\n", *((int64_t *) p->data));
    downsample_i = *((int64_t *) p->data);
    l = g_slice_new0(struct level);
    l->base.downsample = (double) downsample_i;
    l->base.w = data->w / l->base.downsample;
    l->base.h = data->h / l->base.downsample;
    l->downsample_i = downsample_i;

    l->base.tile_w = find_most_common_width(osr, downsample_i);
    l->base.tile_h = find_most_common_height(osr, downsample_i);

    g_ptr_array_add(levels, l);
    p = p->next;
  }

  g_assert(osr->levels == NULL);
  osr->level_count = levels->len;
  osr->levels = (struct _openslide_level **) g_ptr_array_free(levels, false);

  //g_list_free(p);
  g_list_free(downsamples);

  printf("debug: total %d levels\n", osr->level_count);
  return true;
}

static bool load_dir_position(openslide_t *osr, GError **err)
{
  struct zeiss_ops_data *data = osr->data;

  g_autoptr(_openslide_file) f = _openslide_fopen(data->filename, err);
  if (!f)
    return false;

  char buf[512];
  if (!_openslide_fseek(f, sizeof(struct zisraw_seg_hdr), SEEK_SET, err)) {
    g_prefix_error(err, "Couldn't seek to FileHeaderSegment start: ");
    return false;
  }

  size_t len = sizeof(struct zisraw_data_file_hdr);
  if (_openslide_fread(f, buf, len) != len) {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED,
                "Cannot read FileHeader");
    return false;
  }

  struct zisraw_data_file_hdr *hdr = (struct zisraw_data_file_hdr *) buf;
  data->subblk_dir_pos = GINT64_FROM_LE(hdr->subblk_dir_pos);
  data->meta_pos = GINT64_FROM_LE(hdr->meta_pos);
  data->attach_dir_pos = GINT64_FROM_LE(hdr->attach_dir_pos);

  return true;
}

static bool load_czi_structure (openslide_t *osr, GError **err)
{
  //create_db(osr, NULL);
  load_dir_position(osr, err);
  read_subblk_dir(osr, err);
  adjust_coordinate_origin(osr, err);

  return true;
}

static void set_prop(openslide_t *osr, const char *name, const char *value) {
  if (value) {
    g_hash_table_insert(osr->properties,
                        g_strdup(name),
                        g_strdup(value));
  }
}

/* parse XML and set standard openslide properties. Also set width, height in
 * ops_data
 */
static void *parse_xml_set_prop(openslide_t *osr, const char *xml, GError **err)
{
  g_autoptr(xmlDoc) doc = _openslide_xml_parse(xml, err);
  if (doc == NULL) {
    g_printerr("Error: cannot parse XML to xmlDoc\n");
    return NULL;
  }

  // create XPATH context to query the document
  g_autoptr(xmlXPathContext) ctx = _openslide_xml_xpath_create(doc);

  /* part of XML structure:

    ImageDocument
        Metadata
            Experiment
            HardwareSetting
            CustomAttributes
            Information
                User
                Application
                Document
                Image
                    PixelType
                    SizeC
                    SizeX
                    SizeY

                    Dimensions
                        Channels
                            Channel
                            Channel
                        Tracks
                            Track
                            Track
                Instrument
                    Microscopes
                        <Microscope Id="Microscope:1" Name="Axioscan 7">
                    Objectives
                        Objective
                            NominalMagnification  (objective-power)
              Scaling
                  Items
                      <Distance Id="X">  (mpp X)
                          Value  (3.4443237544526617E-07, in meter)
                      <Distance Id="Y">  (mpp Y)
                          Value
  */

  struct zeiss_ops_data *data = osr->data;
  g_autofree char *size_x =
    _openslide_xml_xpath_get_string(ctx,
      "/ImageDocument/Metadata/Information/Image/SizeX/text()");
  data->w = (int32_t) atol(size_x);

  g_autofree char *size_y =
    _openslide_xml_xpath_get_string(ctx,
      "/ImageDocument/Metadata/Information/Image/SizeY/text()");
  data->h = (int32_t) atol(size_y);

  g_autofree char *size_s =
    _openslide_xml_xpath_get_string(ctx,
      "/ImageDocument/Metadata/Information/Image/SizeS/text()");
  data->scene = (int32_t) atol(size_s);

  double d;
  char buf[G_ASCII_DTOSTR_BUF_SIZE];
  // in meter/pixel
  g_autofree char *mpp_x =
    _openslide_xml_xpath_get_string(ctx,
      "/ImageDocument/Metadata/Scaling/Items/Distance[@Id='X']/Value/text()");
  d = _openslide_parse_double(mpp_x);
  g_ascii_dtostr(buf, sizeof(buf), d * 1000000.0);
  // in um/pixel
  set_prop(osr, OPENSLIDE_PROPERTY_NAME_MPP_X, buf);

  g_autofree char *mpp_y =
    _openslide_xml_xpath_get_string(ctx,
      "/ImageDocument/Metadata/Scaling/Items/Distance[@Id='Y']/Value/text()");
  d = _openslide_parse_double(mpp_y);
  g_ascii_dtostr(buf, sizeof(buf), d * 1000000.0);
  set_prop(osr, OPENSLIDE_PROPERTY_NAME_MPP_Y, buf);

  g_autofree char *obj =
    _openslide_xml_xpath_get_string(ctx,
      "/ImageDocument/Metadata/Information/Instrument/Objectives/Objective/NominalMagnification/text()");
  set_prop(osr, OPENSLIDE_PROPERTY_NAME_OBJECTIVE_POWER, obj);

  set_prop(osr, OPENSLIDE_PROPERTY_NAME_VENDOR, "zeiss");
  return NULL;
}

static bool read_czi_meta_xml(openslide_t *osr, GError **err)
{
  struct zeiss_ops_data *data = osr->data;

  g_autoptr(_openslide_file) f = _openslide_fopen(data->filename, err);
  if (!f)
    return false;

  if (!_openslide_fseek(f, data->meta_pos, SEEK_SET, err)) {
    g_prefix_error(err, "Couldn't seek to MetaBlock");
    return false;
  }

  char buf[512];
  size_t len = sizeof(struct zisraw_seg_meta_hdr);
  if (_openslide_fread(f, buf, len) != len) {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED,
                "Cannot read MetaBlock header");
    return false;
  }

  struct zisraw_seg_meta_hdr *hdr = (struct zisraw_seg_meta_hdr *) buf;
  len = (size_t) GINT32_FROM_LE(hdr->xml_size);
  g_autofree char *xml = g_slice_alloc(len + 1);
  if (_openslide_fread(f, xml, len) != len) {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED,
                "Cannot read MetaBlock xml");
    return false;
  }

  xml[len] = '\0';
  parse_xml_set_prop(osr, xml, NULL);

  return true;
}

static bool zeiss_open(openslide_t *osr, const char *filename,
                       struct _openslide_tifflike *t G_GNUC_UNUSED,
                       struct _openslide_hash *quickhash1 G_GNUC_UNUSED,
                       GError **err)
{
  g_autoptr(_openslide_file) f = _openslide_fopen(filename, err);
  if (!f) {
    return false;
  }

  // allocate private data
  struct zeiss_ops_data *data = g_slice_new0(struct zeiss_ops_data);
  data->offset_x = G_MAXINT32;
  data->offset_y = G_MAXINT32;
  data->filename = NULL;
  data->count_tile_width =
    g_hash_table_new_full(g_int64_hash, g_int64_equal,
                          (GDestroyNotify) destroy_int64_key,
                          (GDestroyNotify) destroy_wh_count_hashtable);

  data->count_tile_height =
    g_hash_table_new_full(g_int64_hash, g_int64_equal,
                          (GDestroyNotify) destroy_int64_key,
                          (GDestroyNotify) destroy_wh_count_hashtable);

  // store osr data
  g_assert(osr->data == NULL);
  osr->data = data;
  osr->ops = &zeiss_ops;

  if (osr->levels == NULL) {
    g_mutex_lock(&data->mutex);
    if (!data->filename)
      data->filename = g_strdup(filename);

    if (osr->levels == NULL) {
      load_czi_structure(osr, NULL);
      read_czi_meta_xml(osr, NULL);

      init_levels(osr, err);
      init_range_grids(osr, err);
    }

    g_mutex_unlock(&data->mutex);
  }

  return true;
}

const struct _openslide_format _openslide_format_zeiss = {
  .name = "zeiss",
  .vendor = "zeiss",
  .detect = zeiss_detect,
  .open = zeiss_open,
};
