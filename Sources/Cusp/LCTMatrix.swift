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
    
    public typealias Complex = Numerics.Complex<Scalar>
    
    public var a: Complex
    public var b: Complex
    public var c: Complex
    public var d: Complex
    
    public init(a: Complex, b: Complex, c: Complex, d: Complex) {
        self.a = a
        self.b = b
        self.c = c
        self.d = d
    }
    
}

extension LCTMatrix where Scalar == Double {
    
    public static func identityTransform() -> LCTMatrix {
        return LCTMatrix(a: .one, b: .zero,
                         c: .zero, d: .one)
    }
    
    public static func fourierTransform() -> LCTMatrix {
        return LCTMatrix(a: .zero, b: .one,
                         c: -.one, d: .zero)
    }
    
    public static func inverseFourierTransform() -> LCTMatrix {
        return LCTMatrix(a: .zero, b: -.one,
                         c: .one, d: .zero)
    }
    
    public static func fractionalFourierTransform(angle: Scalar) -> LCTMatrix {
        return LCTMatrix(a: .init(cos(angle)), b: .init(sin(angle)),
                         c: .init(-sin(angle)), d: .init(cos(angle)))
    }
    
}

extension LCTMatrix where Scalar == Double {
    
    public func determinant() -> Complex {
        return a * d - b * c
    }
    
}
