/*
 *  OpenSlide, a library for reading whole slide image files
 *
 *  Copyright (c) 2007-2015 Carnegie Mellon University
 *  Copyright (c) 2011 Google, Inc.
 *  Copyright (c) 2015 Benjamin Gilbert
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

#include <string.h>
#include <config.h>

#include "openslide-private.h"
#include "openslide-decode-jxr.h"

#define BGR24TOARGB32(p)                                                       \
  ((uint32_t)((p)[0]) | ((uint32_t)((p)[1]) << 8) | ((uint32_t)((p)[2]) << 16))

static struct wmp_err_msg {
  ERR id;
  char *msg;
} msgs[] = {
  {WMP_errFail, "WMP_errFail"},
  {WMP_errNotYetImplemented, "WMP_errNotYetImplemented"},
  {WMP_errAbstractMethod, "WMP_errAbstractMethod"},
  {WMP_errOutOfMemory, "WMP_errOutOfMemory"},
  {WMP_errFileIO, "WMP_errFileIO"},
  {WMP_errBufferOverflow, "WMP_errBufferOverflow"},
  {WMP_errInvalidParameter, "WMP_errInvalidParameter"},
  {WMP_errInvalidArgument, "WMP_errInvalidArgument"},
  {WMP_errUnsupportedFormat, "WMP_errUnsupportedFormat"},
  {WMP_errIncorrectCodecVersion, "WMP_errIncorrectCodecVersion"},
  {WMP_errIndexNotFound, "WMP_errIndexNotFound"},
  {WMP_errOutOfSequence, "WMP_errOutOfSequence"},
  {WMP_errNotInitialized, "WMP_errNotInitialized"},
  {WMP_errMustBeMultipleOf16LinesUntilLastCall, "WMP_errMustBeMultipleOf16LinesUntilLastCall"},
  {WMP_errPlanarAlphaBandedEncRequiresTempFile, "WMP_errPlanarAlphaBandedEncRequiresTempFile"},
  {WMP_errAlphaModeCannotBeTranscoded, "WMP_errAlphaModeCannotBeTranscoded"},
  {WMP_errIncorrectCodecSubVersion, "WMP_errIncorrectCodecSubVersion"},
  {0, NULL}
};

static void print_err(ERR err)
{
  if (err >= 0)
    return;

  struct wmp_err_msg *p = &msgs[0];
  fprintf(stderr, "_openslide_jxr_decode_buf error: %ld\n", err);
  while ((p++)->msg) {
    if (p->id == err) {
      fprintf(stderr, "_openslide_jxr_decode_buf error: %s\n", p->msg);
      break;
    }
  }
}

static guint get_bits_per_pixel(const PKPixelFormatGUID *pixel_format)
{
  PKPixelInfo pixel_info;

  pixel_info.pGUIDPixFmt = pixel_format;
  PixelFormatLookup(&pixel_info, LOOKUP_FORWARD);
  return pixel_info.cbitUnit;
}

/* GUID_PKPixelFormat24bppBGR has 24bits per pixel. CAIRO_FORMAT_RGB24 has
 * 32bits, with the upper 8 bits unused
 */
static void convert_24bppbgr_to_cario24bpprgb(struct decoded_jxr *p)
{
  uint32_t new_size = p->w * p->h * 4;
  uint32_t *buf = g_slice_alloc(new_size);
  uint32_t *bp = buf;
  uint32_t i = 0;
  if (p->pixel_size != 24) {
    fprintf(stderr, "Skip convert, pixel size is %d\n", p->pixel_size);
    return;
  }

  while (i < p->size) {
    *bp++ = BGR24TOARGB32(&p->data[i]);
    i += 3;
  }

  g_slice_free1(p->size, p->data);
  p->stride = p->w * 4;
  p->pixel_size = 32;
  p->size = new_size;
  p->data = (uint8_t *) buf;
}

bool _openslide_jxr_decode_buf(void *data, size_t datalen,
                               struct decoded_jxr *dest,
                               GError **unused G_GNUC_UNUSED)
{
  PKFormatConverter *pConverter = NULL;
  PKImageDecode *pDecoder = NULL;
  ERR err = WMP_errSuccess;

  struct WMPStream *pStream = NULL;
  CreateWS_Memory(&pStream, (void *) data, datalen);

  PKPixelFormatGUID fmt;
  // GUID_PKPixelFormat32bppRGBA is not supported by converter
  const PKPixelFormatGUID fmt_out = GUID_PKPixelFormat24bppBGR;
  PKRect rect = {0, 0, 0, 0};

  // IID_PKImageWmpDecode is the only supported decoder PKIID
  Call(PKCodecFactory_CreateCodec(&IID_PKImageWmpDecode, (void **) &pDecoder));
  Call(pDecoder->Initialize(pDecoder, pStream));
  pDecoder->GetSize(pDecoder, &rect.Width, &rect.Height);
  pDecoder->GetPixelFormat(pDecoder, &fmt);

  dest->w = rect.Width;
  dest->h = rect.Height;
  dest->stride = (rect.Width * MAX(get_bits_per_pixel(&fmt),
                                   get_bits_per_pixel(&fmt_out)) + 7) / 8;
  dest->size = dest->stride * dest->h;
  dest->data = g_slice_alloc(dest->size);
  dest->pixel_size = get_bits_per_pixel(&fmt_out);

  //Create color converter
  Call(PKCodecFactory_CreateFormatConverter(&pConverter));
  Call(pConverter->Initialize(pConverter, pDecoder, NULL, fmt_out));
  Call(pConverter->Copy(pConverter, &rect, dest->data, dest->stride));

  pDecoder->Release(&pDecoder);
  pConverter->Release(&pConverter);

  convert_24bppbgr_to_cario24bpprgb(dest);

  //fprintf(stderr, "debug: jxrdecoded w,h = %d,%d\n", dest->w, dest->h);
Cleanup:
  print_err(err);

  return true;
}
