/*
 * Copyright 2015-2016 The OpenSSL Project Authors. All Rights Reserved.
 *
 * Licensed under the OpenSSL license (the "License").  You may not use
 * this file except in compliance with the License.  You can obtain a copy
 * in the file LICENSE in the source distribution or at
 * https://www.openssl.org/source/license.html
 */

#include <stddef.h>

#define POLY1305_BLOCK_SIZE  16
#define POLY1305_DIGEST_SIZE 16
#define POLY1305_KEY_SIZE    32

typedef struct poly1305_context POLY1305;

int GFp_poly1305_init_asm(void *ctx, const unsigned char key[16]);
void GFp_poly1305_blocks(void *ctx, const unsigned char *inp, size_t len,
                     unsigned int padbit);
void GFp_poly1305_emit(void *ctx, unsigned char mac[16],
                   const unsigned int nonce[4]);
