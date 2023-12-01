use subtle::{Choice, CtOption};

use crate::Group;

/// This trait represents an element of a cryptographic group with a large prime-order
/// subgroup and a comparatively-small cofactor.
pub trait Cofactor: Group {
    /// The large prime-order subgroup in which cryptographic operations are performed.
    /// If `Self` implements `PrimeGroup`, then `Self::Subgroup` may be `Self`.
    type Subgroup: Group<Scalar = Self::Scalar> + Into<Self>;

    /// Maps `self` to the prime-order subgroup by multiplying this element by some
    /// `k`-multiple of the cofactor.
    ///
    /// The value `k` does not vary between inputs for a given implementation, but may
    /// vary between different implementations of `CofactorGroup` because some groups have
    /// more efficient methods of clearing the cofactor when `k` is allowed to be
    /// different than `1`.
    ///
    /// If `Self` implements [`PrimeGroup`], this returns `self`.
    fn clear_cofactor(&self) -> Self::Subgroup;

    /// Returns `self` if it is contained in the prime-order subgroup.
    ///
    /// If `Self` implements [`PrimeGroup`], this returns `Some(self)`.
    fn into_subgroup(self) -> CtOption<Self::Subgroup>;

    /// Determines if this element is of small order.
    ///
    /// Returns:
    /// - `true` if `self` is in the torsion subgroup.
    /// - `false` if `self` is not in the torsion subgroup.
    fn is_small_order(&self) -> Choice {
        self.clear_cofactor().is_identity()
    }

    /// Determines if this element is "torsion free", i.e., is contained in the
    /// prime-order subgroup.
    ///
    /// Returns:
    /// - `true` if `self` has trivial torsion and is in the prime-order subgroup.
    /// - `false` if `self` has non-zero torsion component and is not in the prime-order
    ///   subgroup.
    fn is_torsion_free(&self) -> Choice;
}
