//
//  FrFT1D.swift
//  Cusp
//
//  Created by Daniel Clelland on 27/02/2025.
//

import Foundation
import Plinth

extension Matrix where Scalar == Double {
    
    public func frft1D(order: Scalar, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        return ComplexMatrix(real: self).frft1D(order: order, setup: setup)
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    public func frft1D(order: Scalar, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        guard order != 0 else {
            return self
        }
        
        let inputChirp = ComplexMatrix.frftInputChirp(shape: .row(length: shape.length), order: order)
        let outputChirp = ComplexMatrix.frftOutputChirp(shape: .row(length: shape.length), order: order)
        
        let multiplied = self * inputChirp
        let transformed = multiplied.fft1D(setup: setup).fftShifted()
        let result = transformed * outputChirp
        
        return result / Scalar.sqrt(Scalar(shape.count))
    }
    
}
