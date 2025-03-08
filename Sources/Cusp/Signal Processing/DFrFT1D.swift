//
//  DFrFT1D.swift
//  Cusp
//
//  Created by Daniel Clelland on 09/03/2025.
//

import Foundation

extension Matrix where Scalar == Double {
    
    public func dfrft1D(order: Scalar, matrix: ComplexMatrix<Scalar>? = nil) -> ComplexMatrix<Scalar> {
        return ComplexMatrix(real: self).dfrft1D(order: order, matrix: matrix)
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    public func dfrft1D(order: Scalar, matrix: ComplexMatrix<Scalar>? = nil) -> ComplexMatrix<Scalar> {
        let matrix = matrix ?? dfrftMatrix(size: shape.count, order: order)
        return (matrix <*> asColumn()).asRow()
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    public static func dfrftMatrix(size: Int, order: Scalar, approximateOrder: Int = 2) -> ComplexMatrix<Scalar> {
        fatalError("Not implemented yet")
    }
    
}
