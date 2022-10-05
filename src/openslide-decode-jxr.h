/*
 *  OpenSlide, a library for reading whole slide image files
 *
 *  Copyright (c) 2007-2013 Carnegie Mellon University
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
 */

#ifndef OPENSLIDE_OPENSLIDE_DECODE_JXR_H_
#define OPENSLIDE_OPENSLIDE_DECODE_JXR_H_

#include <stdint.h>
#include <glib.h>
#include <JXRGlue.h>

/* JPEG XR support */

typedef struct decoded_jxr {
  uint8_t *data;
  uint32_t w;
  uint32_t h;
  size_t size;
  uint32_t stride;
  uint32_t pixel_size;
} decoded_jxr;

static void decoded_jxr_free(struct decoded_jxr *p) {
  g_slice_free1(p->size, p->data);
  g_slice_free1(sizeof(struct decoded_jxr), p);
}

G_DEFINE_AUTOPTR_CLEANUP_FUNC(decoded_jxr, decoded_jxr_free);

bool _openslide_jxr_decode_buf(void *data, size_t datalen,
                               struct decoded_jxr *dest, GError **err);

#endif
