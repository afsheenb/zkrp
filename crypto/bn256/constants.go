/*
 * Copyright (C) 2019 ING BANK N.V.
 *
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

package bn256

import (
    "github.com/afsheenb/zkrp/util/intconversion"
)

// u is the BN parameter that determines the prime: 1868033³.
var u = intconversion.BigFromBase10("4965661367192848881")

// p is a prime over which we form a basic field: 36u⁴+36u³+24u²+6u+1.
var P = intconversion.BigFromBase10("21888242871839275222246405745257275088696311157297823662689037894645226208583")

// Order is the number of elements in both G₁ and G₂: 36u⁴+36u³+18u²+6u+1.
var Order = intconversion.BigFromBase10("21888242871839275222246405745257275088548364400416034343698204186575808495617")

// xiToPMinus1Over6 is ξ^((p-1)/6) where ξ = i+9.
var xiToPMinus1Over6 = &gfP2{intconversion.BigFromBase10("16469823323077808223889137241176536799009286646108169935659301613961712198316"), intconversion.BigFromBase10("8376118865763821496583973867626364092589906065868298776909617916018768340080")}

// xiToPMinus1Over3 is ξ^((p-1)/3) where ξ = i+9.
var xiToPMinus1Over3 = &gfP2{intconversion.BigFromBase10("10307601595873709700152284273816112264069230130616436755625194854815875713954"), intconversion.BigFromBase10("21575463638280843010398324269430826099269044274347216827212613867836435027261")}

// xiToPMinus1Over2 is ξ^((p-1)/2) where ξ = i+9.
var xiToPMinus1Over2 = &gfP2{intconversion.BigFromBase10("3505843767911556378687030309984248845540243509899259641013678093033130930403"), intconversion.BigFromBase10("2821565182194536844548159561693502659359617185244120367078079554186484126554")}

// xiToPSquaredMinus1Over3 is ξ^((p²-1)/3) where ξ = i+9.
var xiToPSquaredMinus1Over3 = intconversion.BigFromBase10("21888242871839275220042445260109153167277707414472061641714758635765020556616")

// xiTo2PSquaredMinus2Over3 is ξ^((2p²-2)/3) where ξ = i+9 (a cubic root of unity, mod p).
var xiTo2PSquaredMinus2Over3 = intconversion.BigFromBase10("2203960485148121921418603742825762020974279258880205651966")

// xiToPSquaredMinus1Over6 is ξ^((1p²-1)/6) where ξ = i+9 (a cubic root of -1, mod p).
var xiToPSquaredMinus1Over6 = intconversion.BigFromBase10("21888242871839275220042445260109153167277707414472061641714758635765020556617")

// xiTo2PMinus2Over3 is ξ^((2p-2)/3) where ξ = i+9.
var xiTo2PMinus2Over3 = &gfP2{intconversion.BigFromBase10("19937756971775647987995932169929341994314640652964949448313374472400716661030"), intconversion.BigFromBase10("2581911344467009335267311115468803099551665605076196740867805258568234346338")}
