use crate::unknown_order::BigNumber;
use generic_ec::hash_to_curve::Tag;
use libpaillier::{Ciphertext, EncryptionKey, Nonce};
use rand_core::RngCore;

use crate::common::{combine, ProtocolError};
use crate::{EPSILON, L};
pub use crate::common::InvalidProof;

pub struct Data {
    pub q: BigNumber,
    pub key0: EncryptionKey,
    pub c: Ciphertext,
    pub x: BigNumber,
}

pub struct PrivateData {
    pub y: BigNumber,
    pub nonce: Nonce,
}

pub struct Commitment {
    s: BigNumber,
    t: BigNumber,
    a: Ciphertext,
    gamma: BigNumber,
}

pub struct PrivateCommitment {
    alpha: BigNumber,
    mu: BigNumber,
    nu: BigNumber,
    r: Nonce,
}

pub type Challenge = BigNumber;

pub struct Proof {
    z1: BigNumber,
    z2: BigNumber,
    w: BigNumber,
}

pub use crate::common::Aux;

/// Create random commitment
pub fn commit<R: RngCore>(
    aux: &Aux,
    data: &Data,
    pdata: &PrivateData,
    mut rng: R,
) -> Result<(Commitment, PrivateCommitment), ProtocolError> {
    let two_to_l_e = BigNumber::one() << (L + EPSILON);
    let modulo_l = (BigNumber::one() << L) * &aux.rsa_modulo;
    let modulo_l_e = &two_to_l_e * &aux.rsa_modulo;

    let alpha = BigNumber::from_rng(&two_to_l_e, &mut rng);
    let mu = BigNumber::from_rng(&modulo_l, &mut rng);
    let nu = BigNumber::from_rng(&modulo_l_e, &mut rng);

    let (a, r) = data.key0.encrypt(alpha.to_bytes(), None).ok_or(ProtocolError::EncryptionFailed)?;

    let commitment = Commitment {
        s: combine(&aux.s, &pdata.y, &aux.t, &mu, &aux.rsa_modulo),
        t: combine(&aux.s, &alpha, &aux.t, &nu, &aux.rsa_modulo),
        a,
        gamma: &alpha % &data.q,
    };
    let private_commitment = PrivateCommitment {alpha, mu, nu, r};
    Ok((commitment, private_commitment))
}

/// Deterministically compute challenge based on prior known values in protocol
pub fn challenge(
    tag: Tag,
    aux: &Aux,
    data: &Data,
    commitment: &Commitment,
) -> Challenge {
    use sha2::Digest;
    let mut digest = sha2::Sha512::new();

    digest.update(tag.as_bytes());

    digest.update(aux.s.to_bytes());
    digest.update(aux.t.to_bytes());
    digest.update(aux.rsa_modulo.to_bytes());

    digest.update(data.q.to_bytes());
    digest.update(data.key0.to_bytes());
    digest.update(data.c.to_bytes());
    digest.update(data.x.to_bytes());

    digest.update(commitment.s.to_bytes());
    digest.update(commitment.t.to_bytes());
    digest.update(commitment.a.to_bytes());
    digest.update(commitment.gamma.to_bytes());

    // FIXME: hash to bignumber
    BigNumber::from_slice(digest.finalize())
}


pub fn prove(
    data: &Data,
    pdata: &PrivateData,
    pcomm: &PrivateCommitment,
    challenge: &Challenge,
) -> Proof {
    Proof {
        z1: &pcomm.alpha + challenge * &pdata.y,
        z2: &pcomm.nu + challenge * &pcomm.mu,
        w: combine(&pcomm.r, &BigNumber::one(), &pdata.nonce, challenge, data.key0.n()),
    }
}

/// Verify the proof
pub fn verify(
    aux: &Aux,
    data: &Data,
    commitment: &Commitment,
    challenge: &Challenge,
    proof: &Proof,
) -> Result<(), InvalidProof> {
    let one = BigNumber::one();
    fn fail_if(b: bool, msg: InvalidProof) -> Result<(), InvalidProof> {
        if b {
            Ok(())
        } else {
            Err(msg)
        }
    }
    // Three equality checks
    {
        let (lhs, _) = data.key0.encrypt(proof.z1.to_bytes(), Some(proof.w.clone())).ok_or(InvalidProof::EncryptionFailed)?;
        let rhs = combine(&commitment.a, &one, &data.c, challenge, data.key0.nn());
        fail_if(lhs == rhs, InvalidProof::EqualityCheckFailed(1))?;
    }
    {
        let lhs = &proof.z1 % &data.q;
        let rhs = commitment.gamma.modadd(&challenge.modmul(&data.x, &data.q), &data.q);
        fail_if(lhs == rhs, InvalidProof::EqualityCheckFailed(2))?;
    }
    {
        let lhs = combine(&aux.s, &proof.z1, &aux.t, &proof.z2, &aux.rsa_modulo);
        let rhs = combine(&commitment.t, &one, &commitment.s, challenge, &aux.rsa_modulo);
        fail_if(lhs == rhs, InvalidProof::EqualityCheckFailed(3))?;
    }

    Ok(())
}

/// Compute proof for the given data, producing random commitment and
/// deriving determenistic challenge.
///
/// Obtained from the above interactive proof via Fiat-Shamir heuristic.
pub fn compute_proof<R: RngCore>(
    tag: Tag,
    aux: &Aux,
    data: &Data,
    pdata: &PrivateData,
    rng: R,
) -> Result<(Commitment, Challenge, Proof), ProtocolError> {
    let (comm, pcomm) = commit(aux, data, pdata, rng)?;
    let challenge = challenge(tag, aux, data, &comm);
    let proof = prove(data, pdata, &pcomm, &challenge);
    Ok((comm, challenge, proof))
}

#[cfg(test)]
mod test {
    use libpaillier::unknown_order::BigNumber;

    use crate::{EPSILON, L};

    #[test]
    fn passing_test() {
        let private_key0 = libpaillier::DecryptionKey::random().unwrap();
        let key0 = libpaillier::EncryptionKey::from(&private_key0);

        let plaintext = BigNumber::from(28);
        let hiddentext = BigNumber::from(228);
        let modulo = BigNumber::from(100);
        assert_eq!(&plaintext % &modulo, &hiddentext % &modulo);
        let (ciphertext, nonce) = key0.encrypt(plaintext.to_bytes(), None).unwrap();

        let data = super::Data {
            q: modulo,
            key0,
            c: ciphertext,
            x: hiddentext,
        };
        let pdata = super::PrivateData {
            y: plaintext,
            nonce,
        };

        let p = BigNumber::prime(L + EPSILON + 1);
        let q = BigNumber::prime(L + EPSILON + 1);
        let rsa_modulo = p * q;
        let s: BigNumber = 123.into();
        let t: BigNumber = 321.into();
        assert_eq!(s.gcd(&rsa_modulo), 1.into());
        assert_eq!(t.gcd(&rsa_modulo), 1.into());
        let aux = super::Aux { s, t, rsa_modulo };

        let tag = generic_ec::hash_to_curve::Tag::new_unwrap("test".as_bytes());

        let (commitment, challenge, proof) =
            super::compute_proof(tag, &aux, &data, &pdata, rand_core::OsRng::default()).unwrap();
        let r = super::verify(&aux, &data, &commitment, &challenge, &proof);
        match r {
            Ok(()) => (),
            Err(e) => panic!("{:?}", e),
        }
    }
}