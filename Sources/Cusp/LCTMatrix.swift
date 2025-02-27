//
//  LCTMatrix.swift
//  Cusp
//
//  Created by June Russell on 28/02/2025.
//

import Plinth

public struct LCTMatrix<Scalar> {
    
    public var a: Scalar
    public var b: Scalar
    public var c: Scalar
    public var d: Scalar
    
    public init(a: Scalar, b: Scalar, c: Scalar, d: Scalar) {
        self.a = a
        self.b = b
        self.c = c
        self.d = d
    }
    
}

extension LCTMatrix where Scalar == Double {
    
    public func determinant() -> Matrix {
        return a * d - b * c
    }
    
}
