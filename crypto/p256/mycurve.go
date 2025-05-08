/*
 * Original copyright (C) 2019 ING BANK N.V.
 * Rewritten to remove go-ethereum dependency by arbedout
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

package p256

import (
	"crypto/elliptic" // Used for CurveParams which btcec.KoblitzCurve embeds
	"math/big"

	"github.com/btcsuite/btcd/btcec/v2" // WASM-compatible secp256k1 library
)

// MyBitCurve provides an API compatible with the original structure,
// but uses btcec.S256() for its underlying elliptic curve operations
// and parameters.
type MyBitCurve struct {
	// Embed CurveParams to make P, N, B, Gx, Gy, BitSize fields available
	// directly. This matches the structure implicitly provided by the original
	// embedding of go-ethereum's secp256k1.BitCurve.
	// These parameters will be sourced from btcec.S256().CurveParams.
	*elliptic.CurveParams

	// Store the underlying btcec.KoblitzCurve instance to delegate operations to.
	underlyingCurve *btcec.KoblitzCurve
}

// Add returns the sum of (x1,y1) and (x2,y2).
// This version delegates to the btcec library's Add method.
// The original MyBitCurve.Add method called its own addJacobian and affineFromJacobian.
// By delegating to btcec.Add, we use the vetted implementation from the new library.
// The addJacobian and affineFromJacobian methods are preserved below for API compatibility,
// should any external code (not visible here) directly call them.
func (BitCurve *MyBitCurve) Add(x1, y1, x2, y2 *big.Int) (*big.Int, *big.Int) {
	if BitCurve.underlyingCurve == nil {
		// This panic indicates an uninitialized MyBitCurve, which should not occur
		// if S256() is used, as 'theCurve' is initialized in init().
		panic("MyBitCurve not properly initialized: underlyingCurve is nil")
	}
	return BitCurve.underlyingCurve.Add(x1, y1, x2, y2)
}

// addJacobian takes two points in Jacobian coordinates, (x1, y1, z1) and
// (x2, y2, z2) and returns their sum, also in Jacobian form.
// This is the original implementation from mycurve.txt, preserved for API compatibility.
// It relies on BitCurve.P (the field prime), which is now sourced from the
// embedded elliptic.CurveParams (initialized from btcec.S256()).
func (BitCurve *MyBitCurve) addJacobian(x1, y1, z1, x2, y2, z2 *big.Int) (*big.Int, *big.Int, *big.Int) {
	// See http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-add-2007-bl
	// All arithmetic operations are on big.Int and use BitCurve.P for modulo.
	// BitCurve.P is available via the embedded *elliptic.CurveParams.
	z1z1 := new(big.Int).Mul(z1, z1)
	z1z1.Mod(z1z1, BitCurve.P)
	z2z2 := new(big.Int).Mul(z2, z2)
	z2z2.Mod(z2z2, BitCurve.P)

	u1 := new(big.Int).Mul(x1, z2z2)
	u1.Mod(u1, BitCurve.P)
	u2 := new(big.Int).Mul(x2, z1z1)
	u2.Mod(u2, BitCurve.P)
	h := new(big.Int).Sub(u2, u1)
	if h.Sign() == -1 {
		h.Add(h, BitCurve.P)
	}
	i := new(big.Int).Lsh(h, 1)
	i.Mul(i, i)
	j := new(big.Int).Mul(h, i)

	s1 := new(big.Int).Mul(y1, z2)
	s1.Mul(s1, z2z2)
	s1.Mod(s1, BitCurve.P)
	s2 := new(big.Int).Mul(y2, z1)
	s2.Mul(s2, z1z1)
	s2.Mod(s2, BitCurve.P)
	r := new(big.Int).Sub(s2, s1)
	if r.Sign() == -1 {
		r.Add(r, BitCurve.P)
	}
	r.Lsh(r, 1)
	v := new(big.Int).Mul(u1, i)

	x3 := new(big.Int).Set(r)
	x3.Mul(x3, x3)
	x3.Sub(x3, j)
	x3.Sub(x3, v)
	x3.Sub(x3, v)
	x3.Mod(x3, BitCurve.P)

	y3 := new(big.Int).Set(r)
	v.Sub(v, x3)
	y3.Mul(y3, v)
	s1.Mul(s1, j)
	s1.Lsh(s1, 1)
	y3.Sub(y3, s1)
	y3.Mod(y3, BitCurve.P)

	z3 := new(big.Int).Add(z1, z2)
	z3.Mul(z3, z3)
	z3.Sub(z3, z1z1)
	if z3.Sign() == -1 {
		z3.Add(z3, BitCurve.P)
	}
	z3.Sub(z3, z2z2)
	if z3.Sign() == -1 {
		z3.Add(z3, BitCurve.P)
	}
	z3.Mul(z3, h)
	z3.Mod(z3, BitCurve.P)

	return x3, y3, z3
}

// affineFromJacobian converts a point in Jacobian coordinates (x, y, z) to affine form (xOut, yOut).
// This is the original implementation from mycurve.txt, preserved for API compatibility.
// It relies on BitCurve.P (the field prime), which is now sourced from the
// embedded elliptic.CurveParams.
func (BitCurve *MyBitCurve) affineFromJacobian(x, y, z *big.Int) (xOut, yOut *big.Int) {
	if z.Sign() == 0 {
		return new(big.Int), new(big.Int)
	}
	// BitCurve.P is available via the embedded *elliptic.CurveParams.
	zinv := new(big.Int).ModInverse(z, BitCurve.P)
	zinvsq := new(big.Int).Mul(zinv, zinv)

	xOut = new(big.Int).Mul(x, zinvsq)
	xOut.Mod(xOut, BitCurve.P)
	zinvsq.Mul(zinvsq, zinv) // zinvsq now holds zinv^3
	yOut = new(big.Int).Mul(y, zinvsq)
	yOut.Mod(yOut, BitCurve.P)
	return
}

// Double returns 2*(x1,y1).
// This method was previously inherited from the embedded go-ethereum secp256k1.BitCurve.
// It now explicitly delegates to the underlying btcec.KoblitzCurve instance.
func (BitCurve *MyBitCurve) Double(x1, y1 *big.Int) (*big.Int, *big.Int) {
	if BitCurve.underlyingCurve == nil {
		panic("MyBitCurve not properly initialized: underlyingCurve is nil")
	}
	return BitCurve.underlyingCurve.Double(x1, y1)
}

// ScalarMult returns k*(Bx,By).
// This method was previously inherited from the embedded go-ethereum secp256k1.BitCurve.
// It now explicitly delegates to the underlying btcec.KoblitzCurve instance.
// The signature (k []byte) matches the original and btcec's KoblitzCurve.ScalarMult.
func (BitCurve *MyBitCurve) ScalarMult(Bx, By *big.Int, k []byte) (*big.Int, *big.Int) {
	if BitCurve.underlyingCurve == nil {
		panic("MyBitCurve not properly initialized: underlyingCurve is nil")
	}
	return BitCurve.underlyingCurve.ScalarMult(Bx, By, k)
}

// ScalarBaseMult returns k*G, where G is the base point of the curve.
// This method was previously inherited from the embedded go-ethereum secp256k1.BitCurve.
// It now explicitly delegates to the underlying btcec.KoblitzCurve instance.
// The signature (k []byte) matches the original and btcec's KoblitzCurve.ScalarBaseMult.
func (BitCurve *MyBitCurve) ScalarBaseMult(k []byte) (*big.Int, *big.Int) {
	if BitCurve.underlyingCurve == nil {
		panic("MyBitCurve not properly initialized: underlyingCurve is nil")
	}
	return BitCurve.underlyingCurve.ScalarBaseMult(k)
}

// theCurve is the global instance of MyBitCurve, initialized with secp256k1 parameters.
var theCurve = new(MyBitCurve)

func init() {
	// Initialize 'theCurve' using btcec.S256().
	// btcec.S256() returns a *btcec.KoblitzCurve, which embeds *elliptic.CurveParams.
	btcecCurveInstance := btcec.S256()

	// Assign the embedded CurveParams. This makes P, N, B, Gx, Gy, BitSize
	// fields directly accessible on theCurve (e.g., theCurve.P).
	// These parameters are standard for secp256k1 and will match those
	// previously hardcoded or obtained from go-ethereum's library.
	theCurve.CurveParams = btcecCurveInstance.CurveParams

	// Store the btcec.KoblitzCurve instance itself for delegating curve operations.
	theCurve.underlyingCurve = btcecCurveInstance

	// Note: The explicit setting of P, N, B, Gx, Gy, BitSize from the original
	// init() function is no longer needed here, as these are provided by
	// btcecCurveInstance.CurveParams.
}

// S256 returns a MyBitCurve which implements secp256k1 functionalities,
// now backed by the btcec library.
func S256() *MyBitCurve {
	return theCurve
}
