//
//  LCT1D.swift
//  Cusp
//
//  Created by Daniel Clelland on 27/02/2025.
//

import Foundation
import Numerics
import Plinth

extension Matrix where Scalar == Double {
    
    public func lct1D(matrix: LCTMatrix<Scalar>, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        return ComplexMatrix(real: self).lct1D(matrix: matrix, setup: setup)
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    public func lct1D(matrix: LCTMatrix<Scalar>, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        let inputChirp = ComplexMatrix.lctInputChirp(shape: .row(length: shape.length), matrix: matrix)
        let outputChirp = ComplexMatrix.lctOutputChirp(shape: .row(length: shape.length), matrix: matrix)
        
        let multiplied = self * inputChirp
        let transformed = multiplied.fft1D(setup: setup).fftShifted()
        let result = transformed * outputChirp
        
        return result / Scalar.sqrt(Scalar(shape.count))
    }
    
}
