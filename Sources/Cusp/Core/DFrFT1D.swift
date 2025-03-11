//
//  DFrFT1D.swift
//  Cusp
//
//  Created by Daniel Clelland on 09/03/2025.
//

import Foundation
import Plinth

extension Matrix where Scalar == Double {
    
    public func dfrft1D(order: Scalar, derivativeOrder: Int = 2) -> ComplexMatrix<Scalar> {
        return ComplexMatrix(real: self).dfrft1D(order: order, derivativeOrder: derivativeOrder)
    }
    
    public func dfrft1D(matrix: ComplexMatrix<Scalar>) -> ComplexMatrix<Scalar> {
        return ComplexMatrix(real: self).dfrft1D(matrix: matrix)
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    public func dfrft1D(order: Scalar, derivativeOrder: Int = 2) -> ComplexMatrix<Scalar> {
        let matrix = DFrFT<Scalar>(length: shape.count, order: order).matrix(derivativeOrder: derivativeOrder)
        return dfrft1D(matrix: matrix)
    }
    
    public func dfrft1D(matrix: ComplexMatrix<Scalar>) -> ComplexMatrix<Scalar> {
        return (matrix <*> asColumn()).asRow()
    }
    
}
