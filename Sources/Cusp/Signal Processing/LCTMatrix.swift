//
//  LCTMatrix.swift
//  Cusp
//
//  Created by June Russell on 28/02/2025.
//

import Foundation
import Numerics
import Plinth

public struct LCTMatrix<Scalar> where Scalar: Real {
    
    public var a: Scalar
    public var b: Scalar
    public var c: Scalar
    public var d: Scalar
    
    public init(a: Scalar, b: Scalar, c: Scalar, d: Scalar) {
//        let epsilon = 1e-10
//        let determinant = a * d - b * c
//        precondition(abs(determinant - 1.0) < epsilon, "LCTMatrix must be symplectic (ad - bc = 1)")
        
        self.a = a
        self.b = b
        self.c = c
        self.d = d
    }
    
}

extension LCTMatrix where Scalar == Double {
    
    public static let identityTransform = LCTMatrix(
        a: 1.0, b: 0.0,
        c: 0.0, d: 1.0
    )
    
    public static let fourierTransform = LCTMatrix(
        a: 0.0, b: 1.0,
        c: -1.0, d: 0.0
    )
    
    public static let inverseFourierTransform = LCTMatrix(
        a: 0.0, b: -1.0,
        c: 1.0, d: 0.0
    )
    
    public static func fractionalFourierTransform(a: Scalar) -> LCTMatrix {
        let alpha = a * .pi / 2.0
        return LCTMatrix(
            a: cos(alpha), b: sin(alpha),
            c: -sin(alpha), d: cos(alpha)
        )
    }
    
}
