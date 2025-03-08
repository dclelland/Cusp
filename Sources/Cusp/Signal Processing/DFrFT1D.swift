//
//  DFrFT1D.swift
//  Cusp
//
//  Created by Daniel Clelland on 09/03/2025.
//

import Foundation
import Plinth

extension Matrix where Scalar == Double {
    
    public func dfrft1D(order: Scalar, matrix: ComplexMatrix<Scalar>? = nil) -> ComplexMatrix<Scalar> {
        return ComplexMatrix(real: self).dfrft1D(order: order, matrix: matrix)
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    public func dfrft1D(order: Scalar, matrix: ComplexMatrix<Scalar>? = nil) -> ComplexMatrix<Scalar> {
        let matrix = matrix ?? .dfrftMatrix(length: shape.count, order: order)
        return (matrix <*> asColumn()).asRow()
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    public static func dfrftMatrix(length: Int, order: Scalar, approximateOrder: Int = 2) -> ComplexMatrix<Scalar> {
        fatalError("Not implemented yet")
    }
    
}
