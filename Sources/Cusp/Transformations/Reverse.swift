//
//  Reverse.swift
//  Cusp
//
//  Created by Daniel Clelland on 07/03/2025.
//

import Foundation
import Plinth

extension ComplexMatrix where Scalar == Float {
    
    internal func reversedBody() -> ComplexMatrix {
        return ComplexMatrix(real: real.reversedBody(), imaginary: imaginary.reversedBody())
    }
    
}

extension Matrix where Scalar == Float {
    
    internal func reversedBody() -> Matrix {
        let head = elements.first!
        let body = elements.dropFirst()
        return Matrix(shape: shape, elements: [head] + body.reversed())
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    internal func reversedBody() -> ComplexMatrix {
        return ComplexMatrix(real: real.reversedBody(), imaginary: imaginary.reversedBody())
    }
    
}

extension Matrix where Scalar == Double {
    
    internal func reversedBody() -> Matrix {
        let head = elements.first!
        let body = elements.dropFirst()
        return Matrix(shape: shape, elements: [head] + body.reversed())
    }
    
}
