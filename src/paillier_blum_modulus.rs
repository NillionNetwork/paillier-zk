//! ZK-proof of Paillier-Blum modulus. Called Пmod or Rmod in the CGGMP21 paper.
//!
//! ## Description
//! A party P has a modulus `N = pq`, with p and q being Blum primes, and
//! `gcd(N, phi(N)) = 1`. P wants to prove that those equalities about N hold,
//! without disclosing p and q.
//!
//! ## Example
//! 0. Prover P derives two Blum primes and makes a Paillier-Blum modulus
//!     ``` no_run
//!     # use paillier_zk::unknown_order::BigNumber;
//!     fn blum_prime(s: usize) -> BigNumber {
//!         let three = BigNumber::from(3);
//!         loop {
//!             let p = BigNumber::prime(s);
//!             if &p % 4 == three {
//!                 break p;
//!             }
//!         }
//!     }
//!     let p = blum_prime(256);
//!     let q = blum_prime(256);
//!     let n = &p * &q;
//!     // Prover can then make a key from it
//!     let pkey = libpaillier::DecryptionKey::with_primes_unchecked(&p, &q);
//!     ```
//! 1. P computes a non-interactive proof that `n` is a Paillier-Blum modulus:
//!     ``` no_run
//!     use paillier_zk::paillier_blum_modulus as p;
//!     # use generic_ec::hash_to_curve::Tag;
//!     # let (n, p, q) = todo!();
//!     # let shared_state = sha2::Sha256::default();
//!     const TAG: Tag = Tag::new_unwrap("application name".as_bytes());
//!     const SECURITY: usize = 33;
//!
//!     let data = p::Data { n };
//!     let pdata = p::PrivateData { p, q };
//!     let mut rng = rand_core::OsRng::default();
//!
//!     let (commitment, proof) =
//!         p::non_interactive::prove::<{SECURITY}, _, _>(
//!             shared_state,
//!             &data,
//!             &pdata,
//!             &mut rng,
//!         );
//!     ```
//! 2. P sends `data, commitment, proof` to the verifier V
//! 3. V verifies the proof:
//!     ``` no_run
//!     # use generic_ec::hash_to_curve::Tag;
//!     # use paillier_zk::paillier_blum_modulus as p;
//!     # let (data, commitment, proof) = todo!();
//!     # const SECURITY: usize = 33;
//!     # let shared_state = sha2::Sha256::default();
//!     p::non_interactive::verify::<{SECURITY}, _>(
//!         shared_state,
//!         &data,
//!         &commitment,
//!         &proof,
//!     );
//!     ```
//! 4. If the verification succeeded, V can continue communication with P

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

use crate::unknown_order::BigNumber;

/// Reason for failure. If the proof failes, you should only be interested in a
/// reason for debugging purposes
#[derive(Debug, PartialEq, Eq)]
pub enum InvalidProof {
    /// Paillier-Blum modulus is prime
    ModulusIsPrime,
    /// Paillier-Blum modulus
    ModulusIsEven,
    /// Proof's z value in n-th power does not equal commitment value
    IncorrectNthRoot,
    /// Proof's x value in 4-th power does not equal commitment value
    IncorrectFourthRoot,
    /// Couldn't compute modpow
    ModPow,
}

/// Public data that both parties know: the Paillier-Blum modulus
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Data {
    pub n: BigNumber,
}

/// Private data of prover
#[derive(Clone)]
pub struct PrivateData {
    pub p: BigNumber,
    pub q: BigNumber,
}

/// Prover's first message, obtained by [`interactive::commit`]
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Commitment {
    pub w: BigNumber,
}

/// Verifier's challenge to prover. Can be obtained deterministically by
/// [`non_interactive::challenge`] or randomly by [`interactive::challenge`]
///
/// Consists of `M` singular challenges
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct Challenge<const M: usize> {
    pub ys: [BigNumber; M],
}

/// A part of proof. Having enough of those guarantees security
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct ProofPoint {
    pub x: BigNumber,
    pub a: bool,
    pub b: bool,
    pub z: BigNumber,
}

/// The ZK proof. Computed by [`interactive::prove`] or
/// [`non_interactive::prove`]. Consists of M proofs for each challenge
#[derive(Debug, Clone)]
pub struct Proof<const M: usize> {
    pub points: [ProofPoint; M],
}

/// The interactive version of the ZK proof. Should be completed in 3 rounds:
/// prover commits to data, verifier responds with a random challenge, and
/// prover gives proof with commitment and challenge.
pub mod interactive {
    use rand_core::RngCore;

    use crate::common::sqrt::{blum_sqrt, find_residue, non_residue_in};
    use crate::common::BigNumberExt;
    use crate::unknown_order::BigNumber;
    use crate::{Error, ErrorReason};

    use super::{Challenge, Commitment, Data, InvalidProof, PrivateData, Proof, ProofPoint};

    /// Create random commitment
    pub fn commit<R: RngCore>(Data { ref n }: &Data, rng: R) -> Commitment {
        Commitment {
            w: non_residue_in(n, rng),
        }
    }

    /// Compute proof for given data and prior protocol values
    pub fn prove<const M: usize>(
        Data { ref n }: &Data,
        PrivateData { ref p, ref q }: &PrivateData,
        Commitment { ref w }: &Commitment,
        challenge: &Challenge<M>,
    ) -> Result<Proof<M>, Error> {
        let sqrt = |x| blum_sqrt(&x, p, q, n);
        let phi = (p - 1) * (q - 1);
        let n_inverse = n.extended_gcd(&phi).x;
        assert_eq!(n_inverse.modmul(&(n % &phi), &phi), BigNumber::one());
        let points = challenge.ys.clone().map(|y| {
            let z = y.powmod(&n_inverse, n)?;
            let (a, b, y_) = find_residue(&y, w, p, q, n);
            let x = sqrt(sqrt(y_));
            Some(ProofPoint { x, a, b, z })
        });
        if points.iter().any(Option::is_none) {
            return Err(ErrorReason::ModPow.into());
        }
        // we checked that all points are `Some(_)`. Have to do this trick as array::try_map is
        // not yet in stable
        let points = points.map(Option::unwrap);
        Ok(Proof { points })
    }

    /// Verify the proof. If this succeeds, the relation Rmod holds with chance
    /// `1/2^M`
    pub fn verify<const M: usize>(
        data: &Data,
        commitment: &Commitment,
        challenge: &Challenge<M>,
        proof: &Proof<M>,
    ) -> Result<(), InvalidProof> {
        if data.n.is_prime() {
            return Err(InvalidProof::ModulusIsPrime);
        }
        if &data.n % BigNumber::from(2) == BigNumber::zero() {
            return Err(InvalidProof::ModulusIsEven);
        }
        for (point, y) in proof.points.iter().zip(challenge.ys.iter()) {
            if point
                .z
                .powmod(&data.n, &data.n)
                .ok_or(InvalidProof::ModPow)?
                != *y
            {
                return Err(InvalidProof::IncorrectNthRoot);
            }
            let y = y.clone();
            let y = if point.a { &data.n - y } else { y };
            let y = if point.b {
                y.modmul(&commitment.w, &data.n)
            } else {
                y
            };
            if point
                .x
                .powmod(&4.into(), &data.n)
                .ok_or(InvalidProof::ModPow)?
                != y
            {
                return Err(InvalidProof::IncorrectFourthRoot);
            }
        }
        Ok(())
    }

    /// Generate random challenge
    ///
    /// `data` parameter is used to generate challenge in correct range
    pub fn challenge<const M: usize, R: RngCore>(
        Data { ref n }: &Data,
        rng: &mut R,
    ) -> Challenge<M> {
        let ys = [(); M].map(|()| BigNumber::from_rng(n, rng));
        Challenge { ys }
    }
}

/// The non-interactive version of proof. Completed in one round, for example
/// see the documentation of parent module.
pub mod non_interactive {
    use rand_core::RngCore;
    use sha2::{digest::typenum::U32, Digest};

    use super::{Challenge, Commitment, Data, InvalidProof, PrivateData, Proof};
    use crate::unknown_order::BigNumber;
    use crate::Error;

    /// Compute proof for the given data, producing random commitment and
    /// deriving determenistic challenge.
    ///
    /// Obtained from the above interactive proof via Fiat-Shamir heuristic.
    pub fn prove<const M: usize, R: RngCore, D>(
        shared_state: D,
        data: &Data,
        pdata: &PrivateData,
        rng: R,
    ) -> Result<(Commitment, Proof<M>), Error>
    where
        D: Digest<OutputSize = U32> + Clone,
    {
        let commitment = super::interactive::commit(data, rng);
        let challenge = challenge(shared_state, data, &commitment);
        let proof = super::interactive::prove(data, pdata, &commitment, &challenge)?;
        Ok((commitment, proof))
    }

    /// Verify the proof, deriving challenge independently from same data
    pub fn verify<const M: usize, D>(
        shared_state: D,
        data: &Data,
        commitment: &Commitment,
        proof: &Proof<M>,
    ) -> Result<(), InvalidProof>
    where
        D: Digest<OutputSize = U32> + Clone,
    {
        let challenge = challenge(shared_state, data, commitment);
        super::interactive::verify(data, commitment, &challenge, proof)
    }

    /// Deterministically compute challenge based on prior known values in protocol
    pub fn challenge<const M: usize, D>(
        shared_state: D,
        Data { ref n }: &Data,
        commitment: &Commitment,
    ) -> Challenge<M>
    where
        D: Digest<OutputSize = U32> + Clone,
    {
        use rand_core::SeedableRng;
        // since we can't use Default and BigNumber isn't copy, we initialize
        // like this
        let mut ys = [(); M].map(|()| BigNumber::zero());
        for (i, y_ref) in ys.iter_mut().enumerate() {
            let seed = shared_state
                .clone()
                .chain_update(n.to_bytes())
                .chain_update(commitment.w.to_bytes())
                .chain_update((i as u64).to_le_bytes())
                .finalize();
            let mut rng = rand_chacha::ChaCha20Rng::from_seed(seed.into());
            *y_ref = BigNumber::from_rng(n, &mut rng);
        }
        Challenge { ys }
    }
}

#[cfg(test)]
mod test {
    use crate::unknown_order::BigNumber;

    #[test]
    fn passing() {
        let mut rng = rand_dev::DevRng::new();
        let p = blum_prime(256, &mut rng);
        let q = blum_prime(256, &mut rng);
        let n = &p * &q;
        let data = super::Data { n };
        let pdata = super::PrivateData { p, q };
        let shared_state = sha2::Sha256::default();
        let (commitment, proof) = super::non_interactive::prove::<65, _, _>(
            shared_state.clone(),
            &data,
            &pdata,
            &mut rng,
        )
        .unwrap();
        let r = super::non_interactive::verify(shared_state, &data, &commitment, &proof);
        match r {
            Ok(()) => (),
            Err(e) => panic!("{e:?}"),
        }
    }

    #[test]
    fn failing() {
        let mut rng = rand_dev::DevRng::new();
        let p = BigNumber::prime_from_rng(256, &mut rng);
        let q = loop {
            // non blum prime
            let q = BigNumber::prime_from_rng(256, &mut rng);
            if &q % 4 == BigNumber::one() {
                break q;
            }
        };
        let n = &p * &q;
        let data = super::Data { n };
        let pdata = super::PrivateData { p, q };
        let shared_state = sha2::Sha256::default();
        let (commitment, proof) = super::non_interactive::prove::<65, _, _>(
            shared_state.clone(),
            &data,
            &pdata,
            &mut rng,
        )
        .unwrap();
        let r = super::non_interactive::verify(shared_state, &data, &commitment, &proof);
        if r.is_ok() {
            panic!("proof should not pass");
        }
    }

    fn blum_prime<R: rand_core::RngCore>(s: usize, rng: &mut R) -> BigNumber {
        let three = BigNumber::from(3);
        loop {
            let p = BigNumber::prime_from_rng(s, rng);
            if &p % 4 == three {
                break p;
            }
        }
    }
}
