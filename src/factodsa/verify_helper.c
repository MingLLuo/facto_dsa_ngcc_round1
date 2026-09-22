/* Acceptance helper: protocol target from (pk, message) plus the authors'
 * original sig_verify.  The reference source is #included so its static hash
 * routines are reachable without modifying the submitted tree.
 *
 *   verify_helper h      <pkfile> <msgfile>
 *   verify_helper verify <pkfile> <msgfile> <sigfile>
 */
#include <stdio.h>
#include <string.h>
#include "SIG_AlgorithmInstance.c"

DRNG_ctx drng_algorithm;

static int read_file(const char *path, unsigned char *buffer, unsigned long long cap,
                     unsigned long long *len)
{
    FILE *stream = fopen(path, "rb");
    if (!stream) return -1;
    size_t got = fread(buffer, 1, (size_t)cap, stream);
    fclose(stream);
    *len = (unsigned long long)got;
    return 0;
}

int main(int argc, char **argv)
{
    static unsigned char pk[PK_BYTES];
    static unsigned char msg[1 << 20];
    static unsigned char sn[SN_BYTES];
    unsigned long long pk_len = 0, msg_len = 0, sn_len = 0;

    if (argc < 4) {
        fprintf(stderr, "usage: %s h|verify pkfile msgfile [sigfile]\n", argv[0]);
        return 2;
    }
    if (read_file(argv[2], pk, PK_BYTES, &pk_len) != 0 || pk_len != PK_BYTES) {
        fprintf(stderr, "bad public key file\n");
        return 3;
    }
    if (read_file(argv[3], msg, sizeof(msg), &msg_len) != 0) {
        fprintf(stderr, "bad message file\n");
        return 3;
    }
    if (strcmp(argv[1], "h") == 0) {
        unsigned char pkh[FACTO_PKH_BYTES];
        fe_t h[FACTO_M];
        sm3_digest(pk, PK_BYTES, pkh);
        if (xof_field(pkh, msg, msg_len, h) != 0) return 4;
        for (int i = 0; i < FACTO_M; i++)
            printf("%u%c", (unsigned)h[i], i + 1 == FACTO_M ? '\n' : ' ');
        return 0;
    }
    if (strcmp(argv[1], "verify") == 0) {
        if (argc < 5) {
            fprintf(stderr, "signature file required\n");
            return 2;
        }
        if (read_file(argv[4], sn, SN_BYTES, &sn_len) != 0 || sn_len != SN_BYTES) {
            fprintf(stderr, "bad signature length\n");
            return 5;
        }
        int rc = sig_verify(pk, pk_len, sn, sn_len, msg, msg_len);
        printf("%d\n", rc);
        return rc == 0 ? 0 : 1;
    }
    fprintf(stderr, "unknown mode\n");
    return 2;
}
