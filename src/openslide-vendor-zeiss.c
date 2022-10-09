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
#define MAX_CHANNEL       3

#define MERGE_GRAY_TO_ARGB(b, g, r)                                            \
  (0xFF000000 | ((uint32_t)(r) << 16) | ((uint32_t)(g) << 8) | ((uint32_t)(b)))

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
  int64_t att_dir_pos;
};

// SubBlockDirectorySegment
struct zisraw_subblk_dir_hdr {
  struct zisraw_seg_hdr seg_hdr;
  int32_t entry_count;
  char _reserved[124];
  // followed by DirectoryEntryDV list
};

// Metadata segment
struct zisraw_meta_hdr {
  struct zisraw_seg_hdr seg_hdr;
  int32_t xml_size;
  int32_t _attach_size;
  char _reserved[248];
};

// SubBlock segment
struct zisraw_subblk_hdr {
  struct zisraw_seg_hdr seg_hdr;
  int32_t meta_size;
  int32_t attach_size;
  int64_t data_size;
  // followed by DirectoryEntryDV of this subblock
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

// AttachmentEntry - Schema A1
struct __attribute__ ((__packed__)) zisraw_att_entry_a1 {
  char schema[2];
  char _reserved2[10];
  int64_t file_pos;
  int32_t _file_part;
  char guid[16];
  char file_type[8];  // ZIP, ZISRAW, JPG etc.
  char name[80];      // Thumbnail, Label, SlidePreview etc.
};


// Attachment Segment, SID = ZISRAWATTACH
struct __attribute__ ((__packed__)) zisraw_seg_att_hdr {
  struct zisraw_seg_hdr seg_hdr;
  int32_t data_size;
  char _reserved1[12];
  struct zisraw_att_entry_a1 att_entry;
  char _reserved2[112];
  // followed by data
};

// AttachmentDirectory Segment, SID = ZISRAWATTDIR
struct zisraw_att_dir_hdr {
  struct zisraw_seg_hdr seg_hdr;
  int32_t entry_count;
  char _reserved[252];
  // followed by AttachementEntryA1 list
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

// image subblocks from multi-channel with same downsample, x and y store in
// the same subblk array. Process at most 3 channels so they can fit in 32bits
// ARGB output. Leave alpha channel in ARGB untouched.
struct czi_address {
  int64_t downsample;
  int32_t x1, y1;
  uint32_t tw, th;
  struct czi_subblk *subblks[MAX_CHANNEL];
  int32_t len;
  int32_t pixel_type;
  int32_t pixel_bits;
};

struct czi_address_key {
  int64_t x;
  int64_t y;
  int downsample;
};

struct associated_image {
  struct _openslide_associated_image base;
  char *filename;
  int64_t zisraw_offset;
  struct czi_subblk *subblk;
};

static struct associated_image_mapping {
  char *czi_name;
  char *osr_name;
} known_associated_images[] = {
  {"Label", "label"},
  {"SlidePreview", "macro"},
  {NULL, NULL},
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

enum z_pixel_type {
  PT_GRAY8 = 0,
  PT_GRAY16,
  PT_GRAY32FLOAT,
  PT_BGR24,
  PT_BGR48,
  PT_BGR96FLOAT = 8,
  PT_BGRA32,
  PT_GRAY64COMPLEX,
  PT_BGR192COMPLEX,
  PT_GRAY32,
  PT_GRAY64,
};

struct zeiss_ops_data {
  // offset to ZISRAWFILE, one for each file, usually 0. CZI file is like
  // Russian doll, it can embed other CZI files. Non-zero value is the
  // offset to embedded CZI file
  int64_t zisraw_offset;

  char *filename;
  int64_t subblk_dir_pos;
  int64_t meta_pos;
  int64_t att_dir_pos;
  int32_t nsubblk;  // total number of subblocks
  int32_t w;
  int32_t h;
  int32_t pixel_bits;
  int32_t scene;
  int32_t offset_x;
  int32_t offset_y;
  GPtrArray *subblks;
  GHashTable *grids;
  GHashTable *addr_subblks;
  GHashTable *count_levels;
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
  .att_dir_pos = 0,
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

static void destroy_ops_data(struct zeiss_ops_data *data)
{
  g_free(data->filename);
  if (data->count_levels)
    g_hash_table_destroy(data->count_levels);

  if (data->addr_subblks)
    g_hash_table_destroy(data->addr_subblks);

  if (data->grids)
    g_hash_table_destroy(data->grids);

  g_ptr_array_free(data->subblks, TRUE);
  g_slice_free(struct zeiss_ops_data, data);
}

static void destroy(openslide_t *osr) {
  for (int32_t i = 0; i < osr->level_count; i++) {
    destroy_level((struct level *) osr->levels[i]);
  }
  g_free(osr->levels);

  destroy_ops_data(osr->data);
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
  /*
  printf("debug paint_region: ds = %ld, x, y (%ld, %ld), w, h %d/%d\n",
         ds, x, y, w, h);
         */
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
  printf("debug tile: pixeltype%d, comp%d, s%dch%d, ds%ld, x1,y1(%d, %d), w/tw %d/%d, h/th %d/%d, pos %ld\n",
      tile->pixel_type,
      tile->compression,
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
  struct czi_subblk *sb = g_slice_new0(struct czi_subblk);
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

  /*
  if (sb->pixel_type > 0)
    printf("debug read_dir_entry(): %ld %d %d, pixel_type = %d\n",
      sb->downsample_i, sb->x1, sb->y1, sb->pixel_type);
      */

  nread += ndim * 20;
  sb->dir_entry_len = nread;
  g_ptr_array_add(subblks, sb);

  return nread;
}

static bool read_subblk_dir(struct zeiss_ops_data *data, GError **err)
{
  data->subblks = g_ptr_array_new_full(64, (GDestroyNotify) destroy_subblk);

  g_autoptr(_openslide_file) f = _openslide_fopen(data->filename, err);
  if (!f)
    return false;

  char buf[512];
  if (!_openslide_fseek(f, data->zisraw_offset + data->subblk_dir_pos,
                        SEEK_SET, err)) {
    g_prefix_error(err, "Couldn't seek to SubBlock directory");
    return false;
  }

  size_t len;
  len  = sizeof(struct zisraw_subblk_dir_hdr);
  if (_openslide_fread(f, buf, len) != len) {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED,
                "Cannot read FileHeader");
    return false;
  }

  struct zisraw_subblk_dir_hdr *hdr = (struct zisraw_subblk_dir_hdr *) buf;
  if (!g_str_has_prefix(hdr->seg_hdr.sid, "ZISRAWDIRECTORY")) {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED,
                "Not SubBlockDirectory");
    return false;
  }

  data->nsubblk = GINT32_FROM_LE(hdr->entry_count);
  //printf("debug total %d subblocks\n", data->nsubblk);

  len = (size_t) GINT32_FROM_LE(hdr->seg_hdr.allocated_size);
  len -= 128;  // DirectoryEntryDV list starts at offset 128 of segment data
  g_autofree char *buf_dir = g_malloc(len);
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
static bool adjust_coordinate_origin(struct zeiss_ops_data *data,
                                     GError **err G_GNUC_UNUSED)
{
  GPtrArray *subblks = data->subblks;
  struct czi_subblk *b;

  for (guint i = 0; i < subblks->len; i++) {
    b = subblks->pdata[i];
    if (b->x1 < data->offset_x)
      data->offset_x = b->x1;

    if (b->y1 < data->offset_y)
      data->offset_y = b->y1;
  }

  //printf("debug: set offset x,y = %d, %d\n", data->offset_x, data->offset_y);
  for (guint i = 0; i < subblks->len; i++) {
    b = subblks->pdata[i];
    b->x1 -= data->offset_x;
    b->y1 -= data->offset_y;
  }

  return true;
}

static bool czi_uncompressed_read(const char *filename,
                                  int64_t pos, int64_t len, int32_t pixel_type,
                                  int32_t pixel_bits, struct decoded_img *dest,
                                  GError **err)
{
  g_autoptr(_openslide_file) f = _openslide_fopen(filename, err);
  if (!f)
    return false;

  if (!_openslide_fseek(f, pos, SEEK_SET, err)) {
    g_prefix_error(err, "Couldn't seek to jxr pixel data");
    return false;
  }

  dest->size = len;
  dest->data = g_slice_alloc(len);
  //printf("debug czi_uncompressed_read(): len = %ld\n", len);

  if (_openslide_fread(f, dest->data, (size_t) len) != (size_t) len) {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED,
                "Cannot read pixel data");
    g_slice_free1(dest->size, dest->data);
    return false;
  }

  if (pixel_type == PT_BGR48)
    return convert_48bppbgr_to_cario24bpprgb(dest);
  else if (pixel_type == PT_GRAY16)
    return convert_gray16_to_gray8(dest, pixel_bits);
  else if (pixel_type == PT_GRAY8)
    return true;

  return convert_24bppbgr_to_cario24bpprgb(dest);
}





static struct pixel_type_name {
  int id;
  char *fmt;
} pt_names[] = {
	{PT_GRAY8, "GRAY8"},
	{PT_GRAY16, "GRAY16"},
	{PT_BGR24, "BGR24"},
	{PT_GRAY32, "GRAY32"},
	{PT_BGRA32, "BGRA32"},
	{PT_GRAY32FLOAT, "GRAY32FLOAT"},
	{PT_GRAY64, "GRAY64"},
	{PT_GRAY64COMPLEX, "GRAY64COMPLEX"},
	{PT_BGR96FLOAT, "BGR96FLOAT"},
	{PT_BGR192COMPLEX, "BGR192COMPLEX"},
  {-1, NULL}
};


static char *print_pixel_type(int pixel_type)
{
  if (pixel_type < 0)
    return NULL;

  struct pixel_type_name *p = &pt_names[0];
  while ((p++)->fmt) {
    if (p->id == pixel_type)
      return p->fmt;
  }

  return NULL;
}

static bool read_data_from_subblk(const char *filename, int64_t zisraw_offset,
                                  struct czi_subblk *sb,
                                  struct decoded_img *dest, int pixel_bits,
                                  GError **err)
{
  struct zisraw_subblk_hdr *hdr;
  char buf[512];
  size_t len;
  int64_t data_pos, offset_meta, n;

  // work with BGR24, BGR48, GRAY8 and GRAY16
  if (sb->pixel_type != PT_BGR24 && sb->pixel_type != PT_GRAY8 &&
      sb->pixel_type != PT_BGR48 && sb->pixel_type != PT_GRAY16) {
    g_warning("read_data_at_address(): pixel type %d not supported\n",
              sb->pixel_type);
    return false;
  }

  dest->w = sb->tw;
  dest->h = sb->th;
  g_autoptr(_openslide_file) f = _openslide_fopen(filename, err);
  if (!f) {
    g_warning("read_data_from_subblk(): cannot open file %s\n", filename);
    return false;
  }

  if (!_openslide_fseek(f, zisraw_offset + sb->file_pos, SEEK_SET, err)) {
    g_prefix_error(err, "Couldn't seek to SubBlock");
    return false;
  }

  len  = sizeof(struct zisraw_subblk_hdr);
  if (_openslide_fread(f, buf, len) != len) {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED,
                "Cannot read SubBlock header");
    return false;
  }
  hdr = (struct zisraw_subblk_hdr *)buf;
  n = MAX(256 - 16 - sb->dir_entry_len, 0);
  offset_meta = MAX(256, n);
  data_pos = zisraw_offset + sb->file_pos + sizeof(struct zisraw_seg_hdr) +
             offset_meta + GINT32_FROM_LE(hdr->meta_size);
  switch (sb->compression) {
  case COMP_NONE:
    czi_uncompressed_read(filename, data_pos, GINT64_FROM_LE(hdr->data_size),
                          sb->pixel_type, pixel_bits,
                          (struct decoded_img *) dest, NULL);
    break;
  case COMP_JXR:
    _openslide_jxr_read(filename, data_pos, GINT64_FROM_LE(hdr->data_size),
                        pixel_bits, (struct decoded_img *) dest, NULL);
    break;
  case COMP_JPEG:
    g_warning("JPEG is not supported\n");
    return false;
  case COMP_LZW:
    g_warning("LZW is not supported\n");
    return false;
  default:
    g_warning("Unrecognized subblock format\n");
    return false;
  }
  return true;
}

static bool read_data_at_address(const char *filename, int64_t zisraw_offset,
                                 struct czi_address *ca,
                                 struct decoded_img *dest,
                                 GError **err)
{
  struct decoded_img ch[MAX_CHANNEL];
  uint32_t *p;
  size_t gray8_len;
  size_t j = 0;

  // so far only work with BGR24, BGR48, GRAY8 and GRAY16
  if (ca->pixel_type != PT_BGR24 && ca->pixel_type != PT_GRAY8 &&
      ca->pixel_type != PT_BGR48 && ca->pixel_type != PT_GRAY16) {
    g_warning("read_data_at_address(): pixel type %d not supported\n",
              ca->pixel_type);
    return false;
  }

  dest->w = ca->tw;
  dest->h = ca->th;
  if (ca->pixel_type == PT_BGR24) {
    return read_data_from_subblk(filename, zisraw_offset, ca->subblks[0],
                                 dest, ca->pixel_bits, err);
  } else if (ca->pixel_type != PT_GRAY8 && ca->pixel_type != PT_GRAY16) {
    g_warning("read_data_at_address(): pixel type %d not supported\n",
              ca->pixel_type);
    return false;
  }

  for (int i = 0; i < MAX_CHANNEL; i++) {
    if (ca->subblks[i]) {
      ch[i].w = ca->tw;
      ch[i].h = ca->tw;
      // read gray channel
      read_data_from_subblk(filename, zisraw_offset, ca->subblks[i],
                            &(ch[i]), ca->pixel_bits, err);
    } else {
      ch[i].size = ca->tw * ca->th;
      ch[i].data = g_slice_alloc(ch[i].size);
      memset(ch[i].data, 0x7F, ch[i].size);
    }
  }
  // generate a fake RGB from grayscale channels
  dest->size = dest->w * dest->h * 4;
  dest->data = g_slice_alloc(dest->size);
  gray8_len = dest->w * dest->h;
  g_assert(gray8_len == ch[0].size);
  p = (uint32_t *) dest->data;
  while (j < gray8_len) {
    *p++ = MERGE_GRAY_TO_ARGB(ch[0].data[j],
                              ch[1].data[j],
                              ch[2].data[j]);
    j++;
  }

  for (int i = 0; i < MAX_CHANNEL; i++)
    g_slice_free1(ch[i].size, ch[i].data);

  return true;
}

static bool read_tile(openslide_t *osr, cairo_t *cr,
                      struct _openslide_level *level G_GNUC_UNUSED,
                      int64_t tid G_GNUC_UNUSED, void *tile_data,
                      void *arg G_GNUC_UNUSED, GError **err G_GNUC_UNUSED)
{
  struct zeiss_ops_data *data = (struct zeiss_ops_data *) osr->data;
  struct decoded_img dest;
  struct czi_address *ca = tile_data;
  g_autoptr(_openslide_cache_entry) cache_entry = NULL;
  g_autoptr(cairo_surface_t) surface;
  unsigned char *img;
  int stride;

  /*
  fprintf(stderr, "debug read_tile: pos = %ld, ds %ld, x1,y1(%d,%d), tw,th = %d, %d\n",
      sb->file_pos, sb->downsample_i, sb->x1, sb->y1, sb->tw, sb->th);
      */

  // file_pos is unique
  img = (unsigned char *)_openslide_cache_get(osr->cache, 0,
                                               ca->subblks[0]->file_pos, 0,
                                               &cache_entry);
  if (!img) {
    fprintf(stderr, "cache missing\n");
    if (!read_data_at_address(data->filename, data->zisraw_offset,
                              ca, &dest, NULL)) {
      g_warning("read_data_at_address() failed\n");
      return false;
    }
    //g_assert(sb->tw == dest.w);
    //g_assert(sb->th == dest.h);

    // _openslide_cache_entry_unref will free data
    _openslide_cache_put(osr->cache, 0, ca->subblks[0]->file_pos, 0,
                         dest.data, dest.size, &cache_entry);
    img = (unsigned char *)dest.data;
  }

  stride = cairo_format_stride_for_width(CAIRO_FORMAT_RGB24, ca->tw);
  surface = cairo_image_surface_create_for_data(img, CAIRO_FORMAT_RGB24,
                                                ca->tw, ca->th, stride);
  cairo_set_source_surface(cr, surface, 0, 0);
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
  g_free(p);
}

static void count_levels(struct zeiss_ops_data *data, int64_t downsample)
{
  int64_t *k;
  void *unused;

  unused = g_hash_table_lookup(data->count_levels, &downsample);
  if (!unused) {
    k = g_new(int64_t, 1);
    *k = downsample;
    g_hash_table_insert(data->count_levels, k, NULL);
  }
}

static void print_ca(struct czi_address *ca) {
  printf("struct czi_address: ds %ld, x1,y1 (%d, %d) tw,th (%d, %d) len %d, pixel_type %d, pixel_bits %d\n",
      ca->downsample,
      ca->x1,
      ca->y1,
      ca->tw,
      ca->th,
      ca->len,
      ca->pixel_type,
      ca->pixel_bits);
}

static void add_czi_address_to_grid(gpointer k G_GNUC_UNUSED,
                                    gpointer v, gpointer ops_data)
{
  struct _openslide_grid *grid;
  struct czi_address *addr = v;
  struct zeiss_ops_data *data = ops_data;
  /*
  printf("debug add_czi_address_to_grid(): ds %d, pixel_type =%d\n",
     addr->downsample, addr->pixel_type);
  print_ca(addr);
     */

  grid = g_hash_table_lookup(data->grids, &addr->downsample);
  if (!grid)
    return;

  _openslide_grid_range_add_tile(grid,
                                 (double) addr->x1 / addr->downsample,
                                 (double) addr->y1 / addr->downsample,
                                 (double) addr->tw, (double) addr->th, addr);
}

static bool init_range_grids(openslide_t *osr, GError **err G_GNUC_UNUSED)
{
  struct zeiss_ops_data *data = osr->data;
  struct _openslide_grid *grid;
  struct level *l;
  int64_t *k;

  data->grids = g_hash_table_new_full(g_int64_hash, g_int64_equal,
                                      (GDestroyNotify) destroy_int64_key,
                                      (GDestroyNotify) _openslide_grid_destroy);
  for (int i = 0; i < osr->level_count; i++) {
    l = (struct level *) osr->levels[i];
    /*
    fprintf(stderr, "debug: level %d, ds %ld, common tw,th = %ld,%ld\n",
        i, l->downsample_i, l->base.tile_w, l->base.tile_h);
        */
    grid = _openslide_grid_create_range(osr, l->base.tile_w, l->base.tile_h,
                                        read_tile,
                                        NULL);
    k = g_new(int64_t, 1);
    *k = l->downsample_i;
    g_hash_table_insert(data->grids, k, grid);
  }

  g_hash_table_foreach(data->addr_subblks,
                       (GHFunc) add_czi_address_to_grid, data);
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

static bool init_levels(openslide_t *osr, GError **err G_GNUC_UNUSED)
{
  struct zeiss_ops_data *data = osr->data;
  struct czi_subblk *b;
  struct level *l;
  GPtrArray *subblks = data->subblks;
  GPtrArray *levels = g_ptr_array_new();
  int64_t downsample_i;

  for (guint i = 0; i < subblks->len; i++) {
    b = subblks->pdata[i];
    count_levels(data, b->downsample_i);
  }

  GList *downsamples = g_hash_table_get_keys(data->count_levels);
  GList *p = g_list_sort(downsamples, (GCompareFunc) cmp_int64);
  downsamples = p;

  while (p) {
    //printf("debug: downsample_i = %ld\n", *((int64_t *) p->data));
    downsample_i = *((int64_t *) p->data);
    l = g_slice_new0(struct level);
    l->base.downsample = (double) downsample_i;
    l->base.w = data->w / l->base.downsample;
    l->base.h = data->h / l->base.downsample;
    l->downsample_i = downsample_i;
    l->base.tile_w = 256;
    l->base.tile_h = 256;

    g_ptr_array_add(levels, l);
    p = p->next;
  }

  g_assert(osr->levels == NULL);
  osr->level_count = levels->len;
  osr->levels = (struct _openslide_level **) g_ptr_array_free(levels, false);
  g_list_free(downsamples);
  return true;
}

static bool load_dir_position(struct zeiss_ops_data *data, GError **err)
{
  struct zisraw_data_file_hdr *hdr;
  char buf[512];

  g_autoptr(_openslide_file) f = _openslide_fopen(data->filename, err);
  if (!f)
    return false;

  if (!_openslide_fseek(f, data->zisraw_offset + sizeof(struct zisraw_seg_hdr),
                        SEEK_SET, err)) {
    g_prefix_error(err, "Couldn't seek to FileHeaderSegment start: ");
    return false;
  }

  size_t len = sizeof(struct zisraw_data_file_hdr);
  if (_openslide_fread(f, buf, len) != len) {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED,
                "Cannot read FileHeader");
    return false;
  }
  hdr = (struct zisraw_data_file_hdr *) buf;
  data->subblk_dir_pos = GINT64_FROM_LE(hdr->subblk_dir_pos);
  data->meta_pos = GINT64_FROM_LE(hdr->meta_pos);
  data->att_dir_pos = GINT64_FROM_LE(hdr->att_dir_pos);
  return true;
}

static guint czi_addr_key_hash(gconstpointer key)
{
  const struct czi_address_key *k = key;

  return (guint) (k->downsample + 293 * (k->x >> 8) + k->y);
}

static gboolean czi_addr_key_equal(gconstpointer a,
                                               gconstpointer b)
{
  const struct czi_address_key *pa = a;
  const struct czi_address_key *pb = b;

  return (pa->x == pb->x) && (pa->y == pb->y) &&
         (pa->downsample == pb->downsample);
}

static void destroy_addr_key(struct czi_address_key *k) {
  g_slice_free(struct czi_address_key, k);
}

static void destroy_addr(struct czi_address *addr) {
  g_slice_free(struct czi_address, addr);
}

static bool combine_subblks_from_multichannels(struct zeiss_ops_data *data,
                                              GError **err G_GNUC_UNUSED)
{
  struct czi_subblk *b;
  struct czi_address_key key, *key2;
  struct czi_address *addr;
  GPtrArray *subblks = data->subblks;

  data->addr_subblks = g_hash_table_new_full(czi_addr_key_hash,
                                             czi_addr_key_equal,
                                             (GDestroyNotify) destroy_addr_key,
                                             (GDestroyNotify) destroy_addr);
  for (guint i = 0; i < subblks->len; i++) {
    b = subblks->pdata[i];
    if (b->channel >= MAX_CHANNEL)
      continue;

    key.x = b->x1;
    key.y = b->y1;
    key.downsample = b->downsample_i;
    addr = g_hash_table_lookup(data->addr_subblks, &key);
    if (!addr) {
      addr = g_slice_new0(struct czi_address);
      // assuming all channels have same bits, Bgr24 is 8bits, Gray16 usually
      // has 14 bits as described in CZI meta xml
      addr->pixel_bits = data->pixel_bits;
      addr->downsample = b->downsample_i;
      addr->x1 = b->x1;
      addr->y1 = b->y1;
      addr->tw = b->tw;
      addr->th = b->th;
      addr->pixel_type = b->pixel_type;

      key2 = g_slice_new0(struct czi_address_key);
      key2->x = key.x;
      key2->y = key.y;
      key2->downsample = key.downsample;
      g_hash_table_insert(data->addr_subblks, key2, addr);
    }

    if (addr->tw != b->tw || addr->th != b->th ||
        addr->pixel_type != b->pixel_type)
      continue;

    addr->subblks[b->channel] = b;
    addr->len++;
  }
  return true;
}

static void set_prop(openslide_t *osr, const char *name, const char *value)
{
  if (value)
    g_hash_table_insert(osr->properties, g_strdup(name), g_strdup(value));
}

/* parse XML and set standard openslide properties. Also set width, height in
 * ops_data
 */
static void parse_xml_set_prop(openslide_t *osr, const char *xml, GError **err)
{
  struct zeiss_ops_data *data = osr->data;
  double d;
  char buf[G_ASCII_DTOSTR_BUF_SIZE];

  g_autoptr(xmlDoc) doc = _openslide_xml_parse(xml, err);
  if (doc == NULL) {
    g_printerr("Error: cannot parse XML to xmlDoc\n");
    return;
  }

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
                    ComponentBitCount
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
  g_autoptr(xmlXPathContext) ctx = _openslide_xml_xpath_create(doc);

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

  g_autofree char *pixel_bits =
    _openslide_xml_xpath_get_string(ctx,
      "/ImageDocument/Metadata/Information/Image/ComponentBitCount/text()");
  data->pixel_bits = (int32_t) atol(pixel_bits);

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
}

static char *read_czi_meta_xml(struct zeiss_ops_data *data, GError **err)
{
  struct zisraw_meta_hdr *hdr;
  size_t len;
  char buf[512];
  g_autofree char *xml;

  g_autoptr(_openslide_file) f = _openslide_fopen(data->filename, err);
  if (!f)
    return NULL;

  if (!_openslide_fseek(f, data->zisraw_offset + data->meta_pos,
                        SEEK_SET, err)) {
    g_prefix_error(err, "Couldn't seek to MetaBlock");
    return NULL;
  }

  len = sizeof(struct zisraw_meta_hdr);
  if (_openslide_fread(f, buf, len) != len) {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED,
                "Cannot read MetaBlock header");
    return NULL;
  }
  hdr = (struct zisraw_meta_hdr *) buf;

  len = (size_t) GINT32_FROM_LE(hdr->xml_size);
  xml = g_malloc(len + 1);
  if (_openslide_fread(f, xml, len) != len) {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED,
                "Cannot read MetaBlock xml");
    return NULL;
  }
  xml[len] = '\0';
  return g_steal_pointer(&xml);
}

/* find offset to embedded image with @name, such as Label */
static int64_t locate_attachment_by_name(struct zeiss_ops_data *data,
                                         const char *name, GError **err)
{
  struct zisraw_att_dir_hdr *hdr;
  struct zisraw_att_entry_a1 *att;
  int64_t zisraw_offset = 0;
  size_t len;
  int nattch;
  char buf[512];

  g_autoptr(_openslide_file) f = _openslide_fopen(data->filename, err);
  if (!f)
    return false;

  if (!_openslide_fseek(f, data->zisraw_offset + data->att_dir_pos,
                        SEEK_SET, err)) {
    g_prefix_error(err, "Couldn't seek to attachment directory: ");
    return false;
  }

  len = sizeof(struct zisraw_att_dir_hdr);
  if (_openslide_fread(f, buf, len) != len) {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED,
                "Cannot read FileHeader");
    return false;
  }
  hdr = (struct zisraw_att_dir_hdr *) buf;
  nattch = GINT32_FROM_LE(hdr->entry_count);
  //printf("debug total %d attachments\n", nattch);

  len  = sizeof(struct zisraw_att_entry_a1);
  for (int i = 0; i < nattch; i++) {
    if (_openslide_fread(f, buf, len) != len) {
      g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED,
                  "Cannot read attachment directory entry");
      return 0;
    }

    att = (struct zisraw_att_entry_a1 *) buf;
    if (g_strcmp0(att->name, name) == 0) {
      //printf("debug locate_attachment_by_name(): found %s\n", name);
      zisraw_offset = att->file_pos;
      break;
    }
  }

  if (zisraw_offset == 0)
    return 0;

  // + 32 bytes segment header + 256 bytes offset
  return zisraw_offset + 32 + 256;
}

/* @dest is pre-allocated by openslide, size 4 * w * h */
static bool get_associated_image_data(struct _openslide_associated_image *_img,
                                      uint32_t *dest,
                                      GError **err) {
  struct associated_image *img = (struct associated_image *) _img;
  struct decoded_img cbuf;

  //printf("debug get_associated_image_data(): subblk pixel_type = %d\n",
  //       img->subblk->pixel_type);
  read_data_from_subblk(img->filename, img->zisraw_offset, img->subblk,
                        &cbuf, 24, err);
  memcpy(dest, cbuf.data, cbuf.size);
  g_slice_free1(cbuf.size, cbuf.data);
  return true;
}

static void destroy_associated_image(struct _openslide_associated_image *_img) {
  struct associated_image *img = (struct associated_image *) _img;

  g_free(img->filename);
  g_slice_free(struct czi_subblk, img->subblk);
  g_slice_free(struct associated_image, img);
}

static const struct _openslide_associated_image_ops zeiss_associated_ops = {
  .get_argb_data = get_associated_image_data,
  .destroy = destroy_associated_image,
};

static bool _add_associated_image(openslide_t *osr, const char *filename,
                                  const char *name, int64_t zisraw_offset,
                                  struct czi_subblk *sb,
                                  GError **err G_GNUC_UNUSED)
{
  struct associated_image *img = g_slice_new0(struct associated_image);

  img->base.ops = &zeiss_associated_ops;
  img->base.w = sb->tw;
  img->base.h = sb->th;
  img->filename = g_strdup(filename);
  img->zisraw_offset = zisraw_offset;
  img->subblk = g_slice_new(struct czi_subblk);
  memcpy(img->subblk, sb, sizeof(*sb));
  g_hash_table_insert(osr->associated_images, g_strdup(name), img);
  return true;
}

static bool zeiss_add_associated_image(openslide_t *osr, GError **err)
{
  struct zeiss_ops_data *outer_data = (struct zeiss_ops_data *) osr->data;
  struct zeiss_ops_data *data;
  struct associated_image_mapping *map = &known_associated_images[0];
  int64_t zisraw_offset;

  for ( ; map->czi_name; map++) {
    // read the outermost CZI to get offset to ZISRAWFILE
    zisraw_offset = locate_attachment_by_name(outer_data, map->czi_name, err);
    if (zisraw_offset == 0)
      continue;

    data = g_slice_new0(struct zeiss_ops_data);
    data->filename = g_strdup(outer_data->filename);
    data->zisraw_offset = zisraw_offset;

    // knowing offset to ZISRAWFILE, now parse the embeded CZI
    load_dir_position(data, err);
    read_subblk_dir(data, err);
    // expect the embeded CZI file has only one image subblock
    struct czi_subblk *sb = (struct czi_subblk *)data->subblks->pdata[0];
    _add_associated_image(osr, data->filename, map->osr_name,
                          data->zisraw_offset, sb, err);
    destroy_ops_data(data);
  }
  return true;
}

static bool zeiss_open(openslide_t *osr, const char *filename,
                       struct _openslide_tifflike *t G_GNUC_UNUSED,
                       struct _openslide_hash *quickhash1 G_GNUC_UNUSED,
                       GError **err)
{
  /*
  g_autoptr(_openslide_file) f = _openslide_fopen(filename, err);
  if (!f) {
    return false;
  }
  */
  struct zeiss_ops_data *data = g_slice_new0(struct zeiss_ops_data);
  g_assert(osr->data == NULL);
  osr->data = data;
  osr->ops = &zeiss_ops;

  data->zisraw_offset = 0;
  data->offset_x = G_MAXINT32;
  data->offset_y = G_MAXINT32;
  data->filename = g_strdup(filename);
  data->count_levels =
    g_hash_table_new_full(g_int64_hash, g_int64_equal,
                          (GDestroyNotify) destroy_int64_key, NULL);
  load_dir_position(data, NULL);
  read_subblk_dir(data, NULL);
  adjust_coordinate_origin(data, NULL);

  g_autofree char *xml = read_czi_meta_xml(data, NULL);
  parse_xml_set_prop(osr, xml, NULL);
  combine_subblks_from_multichannels(data, NULL);
  init_levels(osr, err);
  init_range_grids(osr, err);
  zeiss_add_associated_image(osr, err);
  return true;
}

const struct _openslide_format _openslide_format_zeiss = {
  .name = "zeiss",
  .vendor = "zeiss",
  .detect = zeiss_detect,
  .open = zeiss_open,
};
