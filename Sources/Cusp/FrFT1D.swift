//
//  FrFT1D.swift
//  Cusp
//
//  Created by Daniel Clelland on 27/02/2025.
//

import Foundation
import Numerics
import Plinth

extension Matrix where Scalar == Double {
    
    public func frft1DMatrix(a: Scalar) -> ComplexMatrix<Scalar> {
        return ComplexMatrix(real: self).frft1DMatrix(a: a)
    }
    
    public func frft1D(a: Scalar, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        return ComplexMatrix(real: self).frft1D(a: a, setup: setup)
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    public func frft1DMatrix(a: Scalar) -> ComplexMatrix {
        guard a != 0 else {
            return self
        }
        
        let kernel = ComplexMatrix.frft1DKernel(shape: shape, a: a)
        let transformed = kernel <*> self.asColumn()
        return transformed.asRow() / Scalar.sqrt(Scalar(shape.count))
    }
    
    public func frft1D(a: Scalar, setup: FFT<Scalar>.Setup? = nil) -> ComplexMatrix<Scalar> {
        // Prepare parameters for the chirp method
        let phi = a * .pi / 2
        
        let N = Scalar(shape.count)
        let normFactor = Scalar.sqrt(N * Scalar.sin(phi).magnitude)
        
        // Create chirp vectors
        let n = Matrix.frftXRamp(shape: .row(length: shape.count))
        
        let cotPhi = 1 / Scalar.tan(phi)
        let chirp1 = .pi * cotPhi * n.square() / N
        let chirpMatrix1 = ComplexMatrix(real: chirp1.cos(), imaginary: chirp1.sin())
        
        let cscPhi = 1 / Scalar.sin(phi)
        let chirp2 = .pi * cscPhi * n.square() / N
        let chirpMatrix2 = ComplexMatrix(real: chirp2.cos(), imaginary: chirp2.sin())
        
        // Step 1: Multiply input by first chirp
        let multiplied = self * chirpMatrix1
        
        // Step 2: Compute FFT
        let transformed = multiplied.fft1D(setup: setup)
        
        // Step 3: Multiply by second chirp
        let result = transformed * (chirpMatrix2)
        
        // Apply normalization
        return result / Scalar.sqrt(Scalar(shape.count))
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    fileprivate static func frft1DKernel(shape: Shape, a: Scalar) -> ComplexMatrix {
        let alpha = a * .pi / 2.0
        let cotAlpha = 1.0 / Scalar.tan(alpha)
        let cscAlpha = 1.0 / Scalar.sin(alpha)
        
        let xRamp: Matrix = .frftXRamp(shape: .square(length: shape.count))
        let yRamp: Matrix = .frftYRamp(shape: .square(length: shape.count))
        
        let quadratic: Matrix = (xRamp.square() + yRamp.square()) * cotAlpha
        let cross = (2.0 * xRamp * yRamp) * cscAlpha
        let phase = (.pi / Scalar(shape.count)) * (quadratic - cross)
        
        return ComplexMatrix(real: phase.cos(), imaginary: phase.sin())
    }
    
}

extension Matrix where Scalar == Double {
        
    fileprivate static func frftXRamp(shape: Shape) -> Matrix {
        let width = shape.columns / 2
        let range = Scalar(-width)...Scalar(width - 1)
        return Matrix.xRamp(shape: shape, range: range)
    }
    
    fileprivate static func frftYRamp(shape: Shape) -> Matrix {
        let height = shape.rows / 2
        let range = Scalar(-height)...Scalar(height - 1)
        return Matrix.yRamp(shape: shape, range: range)
    }
    
}
