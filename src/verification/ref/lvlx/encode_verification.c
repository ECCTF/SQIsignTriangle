#include <verification.h>
#include <string.h>
#include <tutil.h>
#include <fp2.h>
#include <encoded_sizes.h>
#include <assert.h>
#include <torsion_constants.h>

typedef unsigned char byte_t;

// digits

static void
encode_digits(byte_t *enc, const digit_t *x, size_t nbytes)
{
#ifdef TARGET_BIG_ENDIAN
    const size_t ndigits = nbytes / sizeof(digit_t);
    const size_t rem = nbytes % sizeof(digit_t);

    for (size_t i = 0; i < ndigits; i++)
        ((digit_t *)enc)[i] = BSWAP_DIGIT(x[i]);
    if (rem) {
        digit_t ld = BSWAP_DIGIT(x[ndigits]);
        memcpy(enc + ndigits * sizeof(digit_t), (byte_t *)&ld, rem);
    }
#else
    memcpy(enc, (const byte_t *)x, nbytes);
#endif
}

static void
decode_digits(digit_t *x, const byte_t *enc, size_t nbytes, size_t ndigits)
{
    assert(nbytes <= ndigits * sizeof(digit_t));
    memcpy((byte_t *)x, enc, nbytes);
    memset((byte_t *)x + nbytes, 0, ndigits * sizeof(digit_t) - nbytes);

#ifdef TARGET_BIG_ENDIAN
    for (size_t i = 0; i < ndigits; i++)
        x[i] = BSWAP_DIGIT(x[i]);
#endif
}

// fp2_t

static byte_t *
fp2_to_bytes(byte_t *enc, const fp2_t *x)
{
    fp2_encode(enc, x);
    return enc + FP2_ENCODED_BYTES;
}

static const byte_t *
fp2_from_bytes(fp2_t *x, const byte_t *enc)
{
    fp2_decode(x, enc);
    return enc + FP2_ENCODED_BYTES;
}

// curves and points

static byte_t *
proj_to_bytes(byte_t *enc, const fp2_t *x, const fp2_t *z)
{
    assert(!fp2_is_zero(z));
    fp2_t tmp = *z;
    fp2_inv(&tmp);
#ifndef NDEBUG
    {
        fp2_t chk;
        fp2_mul(&chk, z, &tmp);
        fp2_t one;
        fp2_set_one(&one);
        assert(fp2_is_equal(&chk, &one));
    }
#endif
    fp2_mul(&tmp, x, &tmp);
    enc = fp2_to_bytes(enc, &tmp);
    return enc;
}

static const byte_t *
proj_from_bytes(fp2_t *x, fp2_t *z, const byte_t *enc)
{
    enc = fp2_from_bytes(x, enc);
    fp2_set_one(z);
    return enc;
}

static byte_t *
ec_curve_to_bytes(byte_t *enc, const ec_curve_t *curve)
{
    return proj_to_bytes(enc, &curve->A, &curve->C);
}

static const byte_t *
ec_curve_from_bytes(ec_curve_t *curve, const byte_t *enc)
{
    memset(curve, 0, sizeof(*curve));
    return proj_from_bytes(&curve->A, &curve->C, enc);
}

static byte_t *
ec_point_to_bytes(byte_t *enc, const ec_point_t *point)
{
    return proj_to_bytes(enc, &point->x, &point->z);
}

static const byte_t *
ec_point_from_bytes(ec_point_t *point, const byte_t *enc)
{
    return proj_from_bytes(&point->x, &point->z, enc);
}

static byte_t *
ec_basis_to_bytes(byte_t *enc, const ec_basis_t *basis)
{
    enc = ec_point_to_bytes(enc, &basis->P);
    enc = ec_point_to_bytes(enc, &basis->Q);
    enc = ec_point_to_bytes(enc, &basis->PmQ);
    return enc;
}

static const byte_t *
ec_basis_from_bytes(ec_basis_t *basis, const byte_t *enc)
{
    enc = ec_point_from_bytes(&basis->P, enc);
    enc = ec_point_from_bytes(&basis->Q, enc);
    enc = ec_point_from_bytes(&basis->PmQ, enc);
    return enc;
}

//==================================================================================================================================//
static byte_t *
ibz_to_bytes(byte_t *enc, const ibz_t *x, size_t nbytes, bool sgn)
{
#ifndef NDEBUG
    {
        // make sure there is enough space
        ibz_t abs, bnd;
        ibz_init(&bnd);
        ibz_init(&abs);
        ibz_pow(&bnd, &ibz_const_two, 8 * nbytes - sgn);
        ibz_abs(&abs, x);
        assert(ibz_cmp(&abs, &bnd) < 0);
        ibz_finalize(&bnd);
        ibz_finalize(&abs);
    }
#endif
    const size_t digits = (nbytes + sizeof(digit_t) - 1) / sizeof(digit_t);
    digit_t d[digits];
    memset(d, 0, sizeof(d));
    if (ibz_cmp(x, &ibz_const_zero) >= 0) {
        // non-negative, straightforward.
        ibz_to_digits(d, x);
    } else {
        assert(sgn);
        // negative; use two's complement.
        ibz_t tmp;
        ibz_init(&tmp);
        ibz_neg(&tmp, x);
        ibz_sub(&tmp, &tmp, &ibz_const_one);
        ibz_to_digits(d, &tmp);
        for (size_t i = 0; i < digits; ++i)
            d[i] = ~d[i];
#ifndef NDEBUG
        {
            // make sure the result is correct
            ibz_t chk;
            ibz_init(&chk);
            ibz_copy_digit_array(&tmp, d);
            ibz_sub(&tmp, &tmp, x);
            ibz_pow(&chk, &ibz_const_two, 8 * sizeof(d));
            assert(!ibz_cmp(&tmp, &chk));
            ibz_finalize(&chk);
        }
#endif
        ibz_finalize(&tmp);
    }
    encode_digits(enc, d, nbytes);
    return enc + nbytes;
}

static const byte_t *
ibz_from_bytes(ibz_t *x, const byte_t *enc, size_t nbytes, bool sgn)
{
    assert(nbytes > 0);
    const size_t ndigits = (nbytes + sizeof(digit_t) - 1) / sizeof(digit_t);
    assert(ndigits > 0);
    digit_t d[ndigits];
    memset(d, 0, sizeof(d));
    decode_digits(d, enc, nbytes, ndigits);
    if (sgn && enc[nbytes - 1] >> 7) {
        // negative, decode two's complement
        const size_t s = sizeof(digit_t) - 1 - (sizeof(d) - nbytes);
        assert(s < sizeof(digit_t));
        d[ndigits - 1] |= ((digit_t)-1) >> 8 * s << 8 * s;
        for (size_t i = 0; i < ndigits; ++i)
            d[i] = ~d[i];
        ibz_copy_digits(x, d, ndigits);
        ibz_add(x, x, &ibz_const_one);
        ibz_neg(x, x);
    } else {
        // non-negative
        ibz_copy_digits(x, d, ndigits);
    }
    return enc + nbytes;
}
//==================================================================================================================================//

// public API

byte_t *
public_key_to_bytes(byte_t *enc, const public_key_t *pk)
{
#ifndef NDEBUG
    const byte_t *const start = enc;
#endif
    enc = ec_curve_to_bytes(enc, &pk->curve);
    *enc++ = pk->hint_pk;
    assert(enc - start == PUBLICKEY_BYTES);
    return enc;
}

const byte_t *
public_key_from_bytes(public_key_t *pk, const byte_t *enc)
{
#ifndef NDEBUG
    const byte_t *const start = enc;
#endif
    enc = ec_curve_from_bytes(&pk->curve, enc);
    pk->hint_pk = *enc++;
    assert(enc - start == PUBLICKEY_BYTES);
    return enc;
}

void
signature_to_bytes(byte_t *enc, const signature_t *sig)
{
#ifndef NDEBUG
    byte_t *const start = enc;
#endif

    enc = fp2_to_bytes(enc, &sig->E_aux_A);

    *enc++ = sig->backtracking;
    *enc++ = sig->two_resp_length;

    size_t nbytes = (SQIsign_response_length + 9) / 8;
    encode_digits(enc, sig->mat_Bchall_can_to_B_chall[0][0], nbytes);
    enc += nbytes;
    encode_digits(enc, sig->mat_Bchall_can_to_B_chall[0][1], nbytes);
    enc += nbytes;
    encode_digits(enc, sig->mat_Bchall_can_to_B_chall[1][0], nbytes);
    enc += nbytes;
    encode_digits(enc, sig->mat_Bchall_can_to_B_chall[1][1], nbytes);
    enc += nbytes;

    nbytes = SECURITY_BITS / 8;
    encode_digits(enc, sig->chall_coeff, nbytes);
    enc += nbytes;

    *enc++ = sig->hint_aux;
    *enc++ = sig->hint_chall;

    assert(enc - start == SIGNATURE_BYTES);
}

void
signature_from_bytes(signature_t *sig, const byte_t *enc)
{
#ifndef NDEBUG
    const byte_t *const start = enc;
#endif

    enc = fp2_from_bytes(&sig->E_aux_A, enc);

    sig->backtracking = *enc++;
    sig->two_resp_length = *enc++;

    size_t nbytes = (SQIsign_response_length + 9) / 8;
    decode_digits(sig->mat_Bchall_can_to_B_chall[0][0], enc, nbytes, NWORDS_ORDER);
    enc += nbytes;
    decode_digits(sig->mat_Bchall_can_to_B_chall[0][1], enc, nbytes, NWORDS_ORDER);
    enc += nbytes;
    decode_digits(sig->mat_Bchall_can_to_B_chall[1][0], enc, nbytes, NWORDS_ORDER);
    enc += nbytes;
    decode_digits(sig->mat_Bchall_can_to_B_chall[1][1], enc, nbytes, NWORDS_ORDER);
    enc += nbytes;

    nbytes = SECURITY_BITS / 8;
    decode_digits(sig->chall_coeff, enc, nbytes, NWORDS_ORDER);
    enc += nbytes;

    sig->hint_aux = *enc++;
    sig->hint_chall = *enc++;

    assert(enc - start == SIGNATURE_BYTES);
}

//==================================================================================================================================//
void
new_signature_to_bytes(byte_t *enc, const new_signature_t *new_sig)
{
// #ifndef NDEBUG
//     byte_t *const start = enc;
// #endif
//     enc = fp2_to_bytes(enc, &new_sig->E_aux_A);
//     enc = fp2_to_bytes(enc, &new_sig->E_com_A);

//     *enc++ = new_sig->hint_aux;
//     *enc++ = new_sig->c1;
//     *enc++ = new_sig->c2;

//     enc = ibz_to_bytes(enc, &new_sig->q1, NEW_TORSION_2POWER_BYTES, false);
//     enc = ibz_to_bytes(enc, &new_sig->q2, NEW_TORSION_2POWER_BYTES, false);
//     enc = ibz_to_bytes(enc, &new_sig->q3, NEW_TORSION_2POWER_BYTES, true);

//     enc = ibz_to_bytes(enc, &new_sig->mat_Baux_can_to_Baux_image[0][0], TORSION_2POWER_BYTES, false);
//     enc = ibz_to_bytes(enc, &new_sig->mat_Baux_can_to_Baux_image[0][1], TORSION_2POWER_BYTES, false);
//     enc = ibz_to_bytes(enc, &new_sig->mat_Baux_can_to_Baux_image[1][0], TORSION_2POWER_BYTES, false);
//     enc = ibz_to_bytes(enc, &new_sig->mat_Baux_can_to_Baux_image[1][1], TORSION_2POWER_BYTES, false);

//     assert(enc - start == NEW_SIGNATURE_BYTES);
}



void
new_signature_from_bytes(new_signature_t *new_sig, const byte_t *enc)
{





}
//==================================================================================================================================//