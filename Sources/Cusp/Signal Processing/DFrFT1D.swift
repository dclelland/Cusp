//
//  DFrFT1D.swift
//  Cusp
//
//  Created by Daniel Clelland on 09/03/2025.
//

import Foundation

extension Matrix where Scalar == Double {
    
    public func dfrft1D(order: Scalar) -> ComplexMatrix<Scalar> {
        return ComplexMatrix(real: self).dfrft1D(order: order, setup: setup)
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    public func dfrft1D(order: Scalar) -> ComplexMatrix<Scalar> {
        fatalError()
    }
    
}
