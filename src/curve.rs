use alloc::boxed::Box;
use core::fmt;
use core::iter::Sum;
use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use ff::PrimeField;
use rand_core::RngCore;
use subtle::{Choice, ConditionallySelectable, CtOption};

/// A helper trait for types with a group operation.
pub trait GroupOps<Rhs = Self, Output = Self>:
    Add<Rhs, Output = Output> + Sub<Rhs, Output = Output> + AddAssign<Rhs> + SubAssign<Rhs>
{
}

impl<T, Rhs, Output> GroupOps<Rhs, Output> for T where
    T: Add<Rhs, Output = Output> + Sub<Rhs, Output = Output> + AddAssign<Rhs> + SubAssign<Rhs>
{
}

/// A helper trait for references with a group operation.
pub trait GroupOpsOwned<Rhs = Self, Output = Self>: for<'r> GroupOps<&'r Rhs, Output> {}

impl<T, Rhs, Output> GroupOpsOwned<Rhs, Output> for T where T: for<'r> GroupOps<&'r Rhs, Output> {}

/// A helper trait for types implementing group scalar multiplication.
pub trait ScalarMul<Rhs, Output = Self>: Mul<Rhs, Output = Output> + MulAssign<Rhs> {}

impl<T, Rhs, Output> ScalarMul<Rhs, Output> for T where T: Mul<Rhs, Output = Output> + MulAssign<Rhs>
{}

/// A helper trait for references implementing group scalar multiplication.
pub trait ScalarMulOwned<Rhs, Output = Self>: for<'r> ScalarMul<&'r Rhs, Output> {}

impl<T, Rhs, Output> ScalarMulOwned<Rhs, Output> for T where T: for<'r> ScalarMul<&'r Rhs, Output> {}

/// This trait represents an element of a cryptographic group.
pub trait Group:
    Clone + Default + Copy + fmt::Debug + Eq + Sized + Send + Sync + 'static + Neg<Output = Self>
{
    /// Scalars modulo the order of this group's scalar field.
    type Scalar: PrimeField;

    /// Returns an element chosen uniformly at random from the non-identity elements of
    /// this group.
    ///
    /// This function is non-deterministic, and samples from the user-provided RNG.
    fn random(rng: impl RngCore) -> Self;

    /// Returns the additive identity, also known as the "neutral element".
    fn identity() -> Self;

    /// Returns a fixed generator of the prime-order subgroup.
    fn generator() -> Self;

    /// Determines if this point is the identity.
    fn is_identity(&self) -> Choice;
}

/// Affine representation of an elliptic curve point guaranteed to be
/// in the correct prime order subgroup.
pub trait CurveAffine:
    Group
    + Encoding
    + ConditionallySelectable
    + Add<Self, Output = Self::Curve>
    + Sub<Self, Output = Self::Curve>
{
    type Base: ff::Field;
    type Curve: Curve<Affine = Self, Scalar = Self::Scalar, Base = Self::Base>;

    /// Returns whether or not this element is on the curve; should
    /// always be true unless an "unchecked" API was used.
    fn is_on_curve(&self) -> Choice;

    /// Converts this element to its curve representation.
    fn to_curve(&self) -> Self::Curve;

    /// Gets the coordinates of this point.
    ///
    /// Returns None if this is the identity.
    fn coordinates(&self) -> CtOption<Coordinates<Self>>;

    /// Obtains a point given $(x, y)$, failing if it is not on the
    /// curve.
    fn from_coordinates(x: Self::Base, y: Self::Base) -> CtOption<Self>;

    /// Returns the curve constant $a$.
    fn a() -> Self::Base;

    /// Returns the curve constant $b$.
    fn b() -> Self::Base;

    /// Requests a hasher that accepts messages and returns near-uniformly
    /// distributed elements in the group, given domain prefix `domain_prefix`.
    fn hash_to_curve<'a>(domain_prefix: &'a str) -> Box<dyn Fn(&[u8]) -> Self + 'a>;

    /// Apply the curve endomorphism by multiplying the x-coordinate
    /// by an element of multiplicative order 3.
    fn endo(&self) -> Self;
}

/// The affine coordinates of a point on an elliptic curve.
#[cfg(feature = "alloc")]
#[cfg_attr(docsrs, doc(cfg(feature = "alloc")))]
#[derive(Clone, Copy, Debug, Default)]
pub struct Coordinates<C: CurveAffine> {
    pub(crate) x: C::Base,
    pub(crate) y: C::Base,
}

#[cfg(feature = "alloc")]
impl<C: CurveAffine> Coordinates<C> {
    /// Obtains a `Coordinates` value given $(x, y)$, failing if it is not on the curve.
    pub fn new(x: C::Base, y: C::Base) -> CtOption<Self> {
        C::from_coordinates(x, y).map(|_| Coordinates { x, y })
    }
    /// Returns the x-coordinate.
    ///
    /// Equivalent to `Coordinates::u`.
    pub fn x(&self) -> &C::Base {
        &self.x
    }

    /// Returns the y-coordinate.
    ///
    /// Equivalent to `Coordinates::v`.
    pub fn y(&self) -> &C::Base {
        &self.y
    }
}

#[cfg(feature = "alloc")]
impl<C: CurveAffine> ConditionallySelectable for Coordinates<C> {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Coordinates {
            x: C::Base::conditional_select(&a.x, &b.x, choice),
            y: C::Base::conditional_select(&a.y, &b.y, choice),
        }
    }
}

/// Efficient representation of an elliptic curve point guaranteed.
pub trait Curve:
    Group
    + Sum
    + for<'a> Sum<&'a Self>
    + ConditionallySelectable
    + GroupOps<Self::Affine, Self>
    + GroupOps<Self, Self>
    + GroupOpsOwned<Self::Affine>
    + GroupOpsOwned<Self, Self>
    + ScalarMul<<Self as Group>::Scalar>
    + ScalarMulOwned<<Self as Group>::Scalar>
{
    type Base: ff::Field;

    /// Apply the curve endomorphism by multiplying the x-coordinate
    /// by an element of multiplicative order 3.
    fn endo(&self) -> Self;

    /// Returns whether or not this element is on the curve; should
    /// always be true unless an "unchecked" API was used.
    fn is_on_curve(&self) -> Choice;

    /// The affine representation for this elliptic curve.
    type Affine: CurveAffine<Curve = Self, Scalar = Self::Scalar, Base = Self::Base>
        // TODO:
        + Mul<Self::Scalar, Output = Self>
        + for<'r> Mul<&'r Self::Scalar, Output = Self>;

    /// Converts a batch of projective elements into affine elements. This function will
    /// panic if `p.len() != q.len()`.
    fn batch_normalize(p: &[Self], q: &mut [Self::Affine]) {
        assert_eq!(p.len(), q.len());

        for (p, q) in p.iter().zip(q.iter_mut()) {
            *q = p.to_affine();
        }
    }

    /// Gets the coordinates of this point.
    ///
    /// Returns None if this is the identity.
    fn coordinates(&self) -> CtOption<CoordinatesJac<Self>>;
    //    fn jacobian_coordinates(&self) -> (Self::Base, Self::Base, Self::Base);

    /// Obtains a point given $(x, y)$, failing if it is not on the
    /// curve.
    fn from_coordinates(x: Self::Base, y: Self::Base, z: Self::Base) -> CtOption<Self>;
    // fn new_jacobian(x: Self::Base, y: Self::Base, z: Self::Base) -> CtOption<Self>;

    /// Converts this element into its affine representation.
    fn to_affine(&self) -> Self::Affine;

    /// Doubles this element.
    #[must_use]
    fn double(&self) -> Self;
}

/// The affine coordinates of a point on an elliptic curve.
#[cfg(feature = "alloc")]
#[cfg_attr(docsrs, doc(cfg(feature = "alloc")))]
#[derive(Clone, Copy, Debug, Default)]
pub struct CoordinatesJac<C: Curve> {
    pub(crate) x: C::Base,
    pub(crate) y: C::Base,
    pub(crate) z: C::Base,
}

#[cfg(feature = "alloc")]
impl<C: Curve> CoordinatesJac<C> {
    /// Obtains a `Coordinates` value given $(x, y)$, failing if it is not on the curve.
    pub fn new(x: C::Base, y: C::Base, z: C::Base) -> CtOption<Self> {
        C::from_coordinates(x, y, z).map(|_| CoordinatesJac { x, y, z })
    }
    /// Returns the x-coordinate.
    ///
    /// Equivalent to `Coordinates::u`.
    pub fn x(&self) -> &C::Base {
        &self.x
    }

    /// Returns the y-coordinate.
    ///
    /// Equivalent to `Coordinates::v`.
    pub fn y(&self) -> &C::Base {
        &self.y
    }

    /// Returns the z-coordinate.
    ///
    /// Equivalent to `Coordinates::v`.
    pub fn z(&self) -> &C::Base {
        &self.z
    }
}

#[cfg(feature = "alloc")]
impl<C: Curve> ConditionallySelectable for CoordinatesJac<C> {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        CoordinatesJac {
            x: C::Base::conditional_select(&a.x, &b.x, choice),
            y: C::Base::conditional_select(&a.y, &b.y, choice),
            z: C::Base::conditional_select(&a.z, &b.z, choice),
        }
    }
}

pub trait Encoding: Sized {
    /// The encoding of group elements.
    ///
    /// The `Default` implementation is not required to return a valid point encoding. The
    /// bound is present to enable encodings to be constructed generically:
    /// ```
    /// # use group::GroupEncoding;
    /// # use subtle::CtOption;
    /// # struct G;
    /// # impl GroupEncoding for G {
    /// #     type Repr = [u8; 0];
    /// #     fn from_bytes(bytes: &Self::Repr) -> CtOption<Self> { unimplemented!() }
    /// #     fn from_bytes_unchecked(bytes: &Self::Repr) -> CtOption<Self> { unimplemented!() }
    /// #     fn to_bytes(&self) -> Self::Repr { unimplemented!() }
    /// # }
    /// # let buf = &[0u8; 0][..];
    /// let mut encoding = <G as GroupEncoding>::Repr::default();
    /// encoding.as_mut().copy_from_slice(buf);
    /// ```
    ///
    /// It is recommended that the default should be the all-zeroes encoding.
    type Repr: Copy + Default + Send + Sync + 'static + AsRef<[u8]> + AsMut<[u8]>;

    /// Attempts to deserialize a group element from its encoding.
    fn from_bytes(bytes: &Self::Repr) -> CtOption<Self>;

    /// Attempts to deserialize a group element, not checking if the element is valid.
    ///
    /// **This is dangerous to call unless you trust the bytes you are reading; otherwise,
    /// API invariants may be broken.** Please consider using
    /// [`GroupEncoding::from_bytes`] instead.
    fn from_bytes_unchecked(bytes: &Self::Repr) -> CtOption<Self>;

    /// Converts this element into its byte encoding. This may or may not support
    /// encoding the identity.
    // TODO: Figure out how to handle identity encoding generically.
    fn to_bytes(&self) -> Self::Repr;
}

/// Affine representation of a point on an elliptic curve that has a defined uncompressed
/// encoding.
pub trait UncompressedEncoding: Sized {
    type Uncompressed: Default + AsRef<[u8]> + AsMut<[u8]>;

    /// Attempts to deserialize an element from its uncompressed encoding.
    fn from_uncompressed(bytes: &Self::Uncompressed) -> CtOption<Self>;

    /// Attempts to deserialize an uncompressed element, not checking if the element is in
    /// the correct subgroup.
    ///
    /// **This is dangerous to call unless you trust the bytes you are reading; otherwise,
    /// API invariants may be broken.** Please consider using
    /// [`UncompressedEncoding::from_uncompressed`] instead.
    fn from_uncompressed_unchecked(bytes: &Self::Uncompressed) -> CtOption<Self>;

    /// Converts this element into its uncompressed encoding, so long as it's not
    /// the point at infinity.
    fn to_uncompressed(&self) -> Self::Uncompressed;
}
