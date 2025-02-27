//
//  LCT1D.swift
//  Cusp
//
//  Created by Daniel Clelland on 27/02/2025.
//

import Plinth

extension Matrix where Scalar == Double {
    
    public func lct1D(matrix: ComplexMatrix<Scalar>, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        precondition(matrix.shape == .square(length: 2))
        fatalError()
    }
    
}
