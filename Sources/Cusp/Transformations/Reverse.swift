//
//  Reverse.swift
//  Cusp
//
//  Created by Daniel Clelland on 07/03/2025.
//

import Foundation
import Plinth

extension ComplexMatrix where Scalar == Float {
    
    internal func bodyReversed() -> ComplexMatrix {
        return ComplexMatrix(real: real.bodyReversed(), imaginary: imaginary.bodyReversed())
    }
    
}

extension Matrix where Scalar == Float {
    
    internal func bodyReversed() -> Matrix {
        let head = elements.first!
        let body = elements.dropFirst()
        return Matrix(shape: shape, elements: [head] + body.reversed())
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    internal func bodyReversed() -> ComplexMatrix {
        return ComplexMatrix(real: real.bodyReversed(), imaginary: imaginary.bodyReversed())
    }
    
}

extension Matrix where Scalar == Double {
    
    internal func bodyReversed() -> Matrix {
        let head = elements.first!
        let body = elements.dropFirst()
        return Matrix(shape: shape, elements: [head] + body.reversed())
    }
    
}
