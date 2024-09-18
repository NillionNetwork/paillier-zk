//! ZK-proof of paillier multiplication. Called Пmul or Rmul in the CGGMP21 paper.
//!
//! ## Description
//!
//! Both parties P and V have the paillier encryption public key N and
//! the ciphertext-values `X`, `Y` and `C` of
//! the plaintext-values `x`, `y` and `c`, respectively.
//!
//! P wants to show that the plaintext-value `c` of the paillier ciphertext-value `C`
//! is equal to the multiplication of the plaintext-values of `X` and `Y`.
//!
//! P privately holds the plaintext `x`, the `nonce_x` (rho_x in the CGGMP21 paper)
//! and `nonce` (rho in the CGGMP21 paper) such that:
//! 1. `X=enc(x,nonce_x)` is the paillier encryption of `x` with `nonce_x`
//! and public key `N`,
//! 2. `C=(Y^x) * (nonce^N) mod N^2`.
//!
//! Given:
//! - `key`, 'pkey' -  pair of public and private key paillier cryptosystem
//! - `X` - paillier encryption of x with `key`
//! - `Y` - paillier encryption of y with `key`
//! - `C` - paillier encryption of c with `key`
//! - `x` - plaintext - data to obtain proof about
//! - `nonce_x` - rho_x in the paper - data to obtain proof about
//! - `nonce` - rho in the paper - data to obtain proof about
//!
//! Prove:
//! - `decrypt(C) = x * y`
//!
//! Disclosing only: `key`, `X`, `Y`, `C`
//!
//! ## Example
//!
//! ```rust
//! use rug::{Integer, Complete};
//! use paillier_zk::BadExponent;
//! use paillier_zk::{paillier_multiplication as p, IntegerExt};
//! # mod pregenerated {
//! #     use super::*;
//! #     paillier_zk::load_pregenerated_data!(
//! #         prover_decryption_key: fast_paillier::DecryptionKey,
//! #     );
//! # }
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! // Prover and verifier have a shared protocol state
//! let shared_state = "some shared state";
//! let mut rng = rand_core::OsRng;
//! # let mut rng = rand_dev::DevRng::new();
//!
//! // 0. Setup: prover and verifier agree on the level of security:
//!
//! let security = p::SecurityParams {
//!     q: (Integer::ONE << 128_u32).complete(),
//! };
//!
//! // 1. Setup: prover prepares the paillier keys
//!
//! let private_key: fast_paillier::DecryptionKey =
//!     pregenerated::prover_decryption_key();
//! let key = private_key.encryption_key();
//!
//! // 2. Setup: prover has some plaintext `x` and `nonce_x`.
//! // Prover performs a paillier encryption of `x` with `nonce_x` to obtain `X`.
//! // Prover has another `nonce` such that `C=(Y^x) * (nonce^N) mod N^2`.
//!
//! let plaintext_x = Integer::from_rng_pm(key.half_n(), &mut rng);
//! let plaintext_y = Integer::from_rng_pm(key.half_n(), &mut rng);
//! let plaintext_c = (&plaintext_x * plaintext_y.clone()).signed_modulo(key.n());
//! let (X, nonce_x) = key.encrypt_with_random(&mut rng, &plaintext_x)?;
//! let (Y, nonce_y) = key.encrypt_with_random(&mut rng, &plaintext_y).unwrap();
//!
//! // Build the nonce_c
//! let nonce = Integer::gen_invertible(key.n(), &mut rng);
//! let nonce_c = {
//!     let ny_to_px: Integer = nonce_y
//!         .pow_mod_ref(&plaintext_x, key.nn())
//!         .ok_or(BadExponent::undefined()).unwrap()
//!         .into();
//!     key.oadd(&ny_to_px, &nonce).unwrap()
//! }
//!
//! // Encrypt plaintext_c
//! let C = key.encrypt_with(&plaintext_c, &nonce_c).unwrap();
//!
//! // 3. Prover computes a non-interactive proof that the plaintext-value `c`
//! // of the paillier ciphertext-value `C` is equal to
//! // the multiplication of the plaintext-values of `X` and `Y`.
//!
//! let data = p::Data {
//!     key,
//!     ciphertext_x: &X,
//!     ciphertext_y: &Y,
//!     ciphertext_c: &C,
//! };
//! let (commitment, proof) =
//!     p::non_interactive::prove::<sha2::Sha256>(
//!         &shared_state,
//!         data,
//!         p::PrivateData { plaintext_x: &plaintext_x, nonce_x: &nonce_x , nonce : &nonce},
//!         &security,
//!         &mut rng,
//!     )?;
//!
//! // 4. Prover sends this data to verifier
//!
//! # fn send(_: &p::Data, _: &p::Commitment, _: &p::Proof) {  }
//! send(&data, &commitment, &proof);
//!
//! // 5. Verifier receives the data and the proof and verifies it
//!
//! # let recv = || (data, commitment, proof);
//! let (data, commitment, proof) = recv();
//! p::non_interactive::verify::<sha2::Sha256>(
//!     &shared_state,
//!     data,
//!     &commitment,
//!     &security,
//!     &proof,
//! )?;
//! # Ok(()) }
//! ```
//!
//! If the verification succeeded, verifier can continue communication with prover

use fast_paillier::{AnyEncryptionKey, Ciphertext, Nonce};
use rug::Integer;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

pub use crate::common::InvalidProof;

/// Security parameters for proof.
#[derive(Debug, Clone, udigest::Digestable)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SecurityParams {
    /// q in paper. Security parameter for challenge
    #[udigest(as = crate::common::encoding::Integer)]
    pub q: Integer,
}

/// Public data that both parties know
#[derive(Debug, Clone, Copy, udigest::Digestable)]
#[udigest(bound = "")]
pub struct Data<'a> {
    /// N in paper, public key that X was encrypted on
    #[udigest(as = crate::common::encoding::AnyEncryptionKey)]
    pub key: &'a dyn AnyEncryptionKey,
    /// X in paper, x encrypted on N
    #[udigest(as = &crate::common::encoding::Integer)]
    pub ciphertext_x: &'a Ciphertext,
    /// Y in paper, value encrypted on N
    #[udigest(as = &crate::common::encoding::Integer)]
    pub ciphertext_y: &'a Ciphertext,
    /// C in paper, value encrypted on N
    #[udigest(as = &crate::common::encoding::Integer)]
    pub ciphertext_c: &'a Ciphertext,
}

/// Private data of prover
#[derive(Clone, Copy)]
pub struct PrivateData<'a> {
    /// x in paper, and plaintext of X
    pub plaintext_x: &'a Integer,
    /// rho_x in paper, nonce_x in encryption x -> X
    pub nonce_x: &'a Nonce,
    /// rho in paper, nonce in relation C<->Y
    pub nonce: &'a Nonce,
}

/// Prover's first message, obtained by [`interactive::commit`]
#[derive(Debug, Clone, udigest::Digestable)]
#[udigest(bound = "")]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize), serde(bound = ""))]
pub struct Commitment {
    #[udigest(as = crate::common::encoding::Integer)]
    pub a: Integer,
    #[udigest(as = crate::common::encoding::Integer)]
    pub b: Integer,
}

/// Prover's data accompanying the commitment. Kept as state between rounds in
/// the interactive protocol.
#[derive(Clone)]
pub struct PrivateCommitment {
    pub alpha: Integer,
    pub r: Integer,
    pub s: Integer,
}

/// Verifier's challenge to prover. Can be obtained deterministically by
/// [`non_interactive::challenge`] or randomly by [`interactive::challenge`]
pub type Challenge = Integer;

/// The ZK proof. Computed by [`interactive::prove`] or
/// [`non_interactive::prove`]
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Proof {
    pub z: Integer,
    pub u: Integer,
    pub v: Integer,
}

/// The interactive version of the ZK proof. Should be completed in 3 rounds:
/// prover commits to data, verifier responds with a random challenge, and
/// prover gives proof with commitment and challenge.
pub mod interactive {
    use rand_core::RngCore;
    use rug::{Complete, Integer};

    use crate::common::{fail_if_ne, IntegerExt, InvalidProofReason};
    use crate::{BadExponent, Error, InvalidProof};

    use super::{
        Challenge, Commitment, Data, PrivateCommitment, PrivateData, Proof, SecurityParams,
    };

    /// Create random commitment
    pub fn commit<R: RngCore>(
        data: Data,
        mut rng: R,
    ) -> Result<(Commitment, PrivateCommitment), Error> {
        let alpha = Integer::from_rng_pm(data.key.half_n(), &mut rng);
        let r = Integer::gen_invertible(data.key.n(), &mut rng);
        let s = Integer::gen_invertible(data.key.n(), &mut rng);

        let a = {
            let alpha_at_y = data.key.omul(&alpha, data.ciphertext_y)?;
            let r_to_key: Integer = r
                .pow_mod_ref(data.key.n(), data.key.nn())
                .ok_or(BadExponent::undefined())?
                .into();
            (alpha_at_y * r_to_key).modulo(data.key.nn())
        };
        let b = data.key.encrypt_with(&alpha, &s)?;

        let commitment = Commitment { a, b };
        let private_commitment = PrivateCommitment { alpha, r, s };
        Ok((commitment, private_commitment))
    }

    /// Compute proof for given data and prior protocol values
    pub fn prove(
        data: Data,
        pdata: PrivateData,
        pcomm: &PrivateCommitment,
        challenge: &Challenge,
    ) -> Result<Proof, Error> {
        let nonce_to_challenge_mod_n: Integer = pdata
            .nonce
            .pow_mod_ref(challenge, data.key.n())
            .ok_or(BadExponent::undefined())?
            .into();
        let nonce_x_to_challenge_mod_n: Integer = pdata
            .nonce_x
            .pow_mod_ref(challenge, data.key.n())
            .ok_or(BadExponent::undefined())?
            .into();
        Ok(Proof {
            z: (&pcomm.alpha + (challenge * pdata.plaintext_x)).complete(),
            u: (&pcomm.r * nonce_to_challenge_mod_n).modulo(data.key.n()),
            v: (&pcomm.s * nonce_x_to_challenge_mod_n).modulo(data.key.n()),
        })
    }

    /// Verify the proof
    pub fn verify(
        data: Data,
        commitment: &Commitment,
        challenge: &Challenge,
        proof: &Proof,
    ) -> Result<(), InvalidProof> {
        {
            // lhs = (1 + N)^z * v^N mod N^2 = Enc(z, v)
            let lhs = data
                .key
                .encrypt_with(&proof.z.signed_modulo(data.key.n()), &proof.v)
                .map_err(|_| InvalidProofReason::PaillierEnc)?;
            let rhs = {
                let e_at_x = data
                    .key
                    .omul(challenge, data.ciphertext_x)
                    .map_err(|_| InvalidProofReason::PaillierOp)?;
                data.key
                    .oadd(&commitment.b, &e_at_x)
                    .map_err(|_| InvalidProofReason::PaillierOp)?
            };
            fail_if_ne(InvalidProofReason::EqualityCheck(1), lhs, rhs)?;
        }
        {
            let lhs = {
                let z_at_y = data
                    .key
                    .omul(&proof.z, data.ciphertext_y)
                    .map_err(|_| InvalidProofReason::PaillierOp)?;
                let u_to_key: Integer = proof
                    .u
                    .pow_mod_ref(data.key.n(), data.key.nn())
                    .ok_or(BadExponent::undefined())?
                    .into();
                (z_at_y * u_to_key).modulo(data.key.nn())
            };
            let rhs = {
                let e_at_c = data
                    .key
                    .omul(challenge, data.ciphertext_c)
                    .map_err(|_| InvalidProofReason::PaillierOp)?;
                data.key
                    .oadd(&commitment.a, &e_at_c)
                    .map_err(|_| InvalidProofReason::PaillierOp)?
            };
            fail_if_ne(InvalidProofReason::EqualityCheck(2), lhs, rhs)?;
        }
        Ok(())
    }

    /// Generate random challenge
    ///
    /// `data` parameter is used to generate challenge in correct range
    pub fn challenge<R>(security: &SecurityParams, rng: &mut R) -> Integer
    where
        R: RngCore,
    {
        Integer::from_rng_pm(&security.q, rng)
    }
}

/// The non-interactive version of proof. Completed in one round, for example
/// see the documentation of parent module.
pub mod non_interactive {
    use digest::Digest;
    use rand_core::RngCore;

    use crate::{Error, InvalidProof};

    use super::{Challenge, Commitment, Data, PrivateData, Proof, SecurityParams};

    /// Compute proof for the given data, producing random commitment and
    /// deriving determenistic challenge.
    ///
    /// Obtained from the above interactive proof via Fiat-Shamir heuristic.
    pub fn prove<D: Digest>(
        shared_state: &impl udigest::Digestable,
        data: Data,
        pdata: PrivateData,
        security: &SecurityParams,
        rng: &mut impl RngCore,
    ) -> Result<(Commitment, Proof), Error> {
        let (comm, pcomm) = super::interactive::commit(data, rng)?;
        let challenge = challenge::<D>(shared_state, data, &comm, security);
        let proof = super::interactive::prove(data, pdata, &pcomm, &challenge)?;
        Ok((comm, proof))
    }

    /// Verify the proof, deriving challenge independently from same data
    pub fn verify<D: Digest>(
        shared_state: &impl udigest::Digestable,
        data: Data,
        commitment: &Commitment,
        security: &SecurityParams,
        proof: &Proof,
    ) -> Result<(), InvalidProof> {
        let challenge = challenge::<D>(shared_state, data, commitment, security);
        super::interactive::verify(data, commitment, &challenge, proof)
    }

    /// Internal function for deriving challenge from protocol values
    /// deterministically
    pub fn challenge<D: Digest>(
        shared_state: &impl udigest::Digestable,
        data: Data,
        commitment: &Commitment,
        security: &SecurityParams,
    ) -> Challenge {
        let tag = "paillier-zk.paillier_multiplication.ni-challenge";
        let seed = udigest::inline_struct!(tag {
            shared_state,
            security,
            data,
            commitment,
        });
        let mut rng = rand_hash::HashRng::<D, _>::from_seed(seed);
        super::interactive::challenge(security, &mut rng)
    }
}

#[cfg(test)]
mod test {
    use rug::{Complete, Integer};
    use sha2::Digest;

    use crate::common::{IntegerExt, InvalidProofReason};
    use crate::BadExponent;

    fn run_with<D: Digest>(
        rng: &mut impl rand_core::CryptoRngCore,
        security: super::SecurityParams,
        data: super::Data,
        pdata: super::PrivateData,
    ) -> Result<(), crate::common::InvalidProof> {
        let shared_state = "shared state";
        let (commitment, proof) =
            super::non_interactive::prove::<D>(&shared_state, data, pdata, &security, rng).unwrap();
        super::non_interactive::verify::<D>(&shared_state, data, &commitment, &security, &proof)
    }

    #[test]
    fn passing() {
        let mut rng = rand_dev::DevRng::new();
        let security = super::SecurityParams {
            q: (Integer::ONE << 128_u32).complete(),
        };
        let private_key = crate::common::test::random_key(&mut rng).unwrap();
        let key = private_key.encryption_key();

        let plaintext_x = Integer::from_rng_pm(key.half_n(), &mut rng);
        let plaintext_y = Integer::from_rng_pm(key.half_n(), &mut rng);
        let plaintext_c = (&plaintext_x * plaintext_y.clone()).signed_modulo(key.n());
        let (ciphertext_x, nonce_x) = key.encrypt_with_random(&mut rng, &plaintext_x).unwrap();
        let (ciphertext_y, nonce_y) = key.encrypt_with_random(&mut rng, &plaintext_y).unwrap();

        // Build nonce_c
        let nonce = Integer::gen_invertible(key.n(), &mut rng);
        let nonce_c = {
            let ny_to_px: Integer = nonce_y
                .pow_mod_ref(&plaintext_x, key.nn())
                .ok_or(BadExponent::undefined())
                .unwrap()
                .into();
            key.oadd(&ny_to_px, &nonce).unwrap()
        };
        // Encrypt c
        let ciphertext_c = key.encrypt_with(&plaintext_c, &nonce_c).unwrap();

        // Create input data
        let data = super::Data {
            key,
            ciphertext_x: &ciphertext_x,
            ciphertext_y: &ciphertext_y,
            ciphertext_c: &ciphertext_c,
        };
        let pdata = super::PrivateData {
            plaintext_x: &plaintext_x,
            nonce_x: &nonce_x,
            nonce: &nonce,
        };

        let r = run_with::<sha2::Sha256>(&mut rng, security, data, pdata);
        match r {
            Ok(()) => (),
            Err(e) => panic!("{e:?}"),
        }
    }
    #[test]
    fn failing_eq_check_1() {
        // Scenario where prover P does not know plaintext x
        let mut rng = rand_dev::DevRng::new();
        let security = super::SecurityParams {
            q: (Integer::ONE << 128_u32).complete(),
        };
        let private_key = crate::common::test::random_key(&mut rng).unwrap();
        let key = private_key.encryption_key();

        // Generate random plaintext for x
        let plaintext_x_random = Integer::from_rng_pm(key.half_n(), &mut rng);

        // Generate correctly C
        let unkown_plaintext_x = Integer::from_rng_pm(key.half_n(), &mut rng);
        let plaintext_y = Integer::from_rng_pm(key.half_n(), &mut rng);
        let plaintext_c = (&unkown_plaintext_x * plaintext_y.clone()).signed_modulo(key.n());
        let (ciphertext_x, nonce_x) = key
            .encrypt_with_random(&mut rng, &unkown_plaintext_x)
            .unwrap();
        let (ciphertext_y, nonce_y) = key.encrypt_with_random(&mut rng, &plaintext_y).unwrap();
        // Build nonce_c
        let nonce = Integer::gen_invertible(key.n(), &mut rng);
        let nonce_c = {
            let ny_to_upx: Integer = nonce_y
                .pow_mod_ref(&unkown_plaintext_x, key.nn())
                .ok_or(BadExponent::undefined())
                .unwrap()
                .into();
            key.oadd(&ny_to_upx, &nonce).unwrap()
        };
        // Encrypt c
        let ciphertext_c = key.encrypt_with(&plaintext_c, &nonce_c).unwrap();

        // Create input data
        let data = super::Data {
            key,
            ciphertext_x: &ciphertext_x,
            ciphertext_y: &ciphertext_y,
            ciphertext_c: &ciphertext_c,
        };
        let pdata = super::PrivateData {
            plaintext_x: &plaintext_x_random,
            nonce_x: &nonce_x,
            nonce: &nonce,
        };

        let r = run_with::<sha2::Sha256>(&mut rng, security, data, pdata)
            .expect_err("proof should not pass");
        match r.reason() {
            InvalidProofReason::EqualityCheck(1) => (),
            e => panic!("proof should not fail with {e:?}"),
        }
    }
    #[test]
    fn failing_eq_check_2() {
        // Scenario where prover P does not know nonce
        let mut rng = rand_dev::DevRng::new();
        let security = super::SecurityParams {
            q: (Integer::ONE << 128_u32).complete(),
        };
        let private_key = crate::common::test::random_key(&mut rng).unwrap();
        let key = private_key.encryption_key();

        // Generate random elements
        let plaintext_x = Integer::from_rng_pm(key.half_n(), &mut rng);
        let plaintext_y = Integer::from_rng_pm(key.half_n(), &mut rng);
        let plaintext_c = Integer::from_rng_pm(key.half_n(), &mut rng);
        let (ciphertext_x, nonce_x) = key.encrypt_with_random(&mut rng, &plaintext_x).unwrap();
        let (ciphertext_y, _) = key.encrypt_with_random(&mut rng, &plaintext_y).unwrap();
        let (ciphertext_c, nonce_c) = key.encrypt_with_random(&mut rng, &plaintext_c).unwrap();

        // Create input data
        let data = super::Data {
            key,
            ciphertext_x: &ciphertext_x,
            ciphertext_y: &ciphertext_y,
            ciphertext_c: &ciphertext_c,
        };
        let pdata = super::PrivateData {
            plaintext_x: &plaintext_x,
            nonce_x: &nonce_x,
            nonce: &nonce_c,
        };

        let r = run_with::<sha2::Sha256>(&mut rng, security, data, pdata)
            .expect_err("proof should not pass");
        match r.reason() {
            InvalidProofReason::EqualityCheck(2) => (),
            e => panic!("proof should not fail with {e:?}"),
        }
    }
}
