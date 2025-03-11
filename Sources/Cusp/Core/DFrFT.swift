//
//  DFrFT.swift
//  Cusp
//
//  Created by Daniel Clelland on 11/03/2025.
//

import Foundation
import Numerics
import Plinth

public struct DFrFT<Scalar> where Scalar: Real {
    
    public var length: Int
    public var order: Scalar
    
    public init(length: Int, order: Scalar) {
        self.length = length
        self.order = order
    }
    
}

extension DFrFT where Scalar == Double {
    
    public func matrix(derivativeOrder: Int = 2) -> ComplexMatrix<Scalar> {
        let eigenvectors = eigenvectors(derivativeOrder: derivativeOrder)
        let eigenvalues = eigenvalues()
        return eigenvectors <*> .diagonal(vector: eigenvalues.elements) <*> eigenvectors.conjugate().transposed()
    }
    
    public func eigenvectors(derivativeOrder: Int = 2) -> ComplexMatrix<Scalar> {
        let hamiltonian = hamiltonian(derivativeOrder: derivativeOrder)
        let decomposition = decomposition()
        let decomposed: Matrix = decomposition <*> hamiltonian <*> decomposition.transposed()
        
        let cosines: Matrix = decomposed[0..<(length / 2 + 1), 0..<(length / 2 + 1)]
        let sines: Matrix = decomposed[(length / 2 + 1)..<length, (length / 2 + 1)..<length]
        
        let cosineEigenvectors = try! cosines.eigendecomposition(computing: .rightEigenvectors).sorted(.realAscending).rightEigenvectors.real
        let sineEigenvectors = try! sines.eigendecomposition(computing: .rightEigenvectors).sorted(.realAscending).rightEigenvectors.real
        
        let transformedCosineEigenvectors = (decomposition <*> cosineEigenvectors.padded(bottom: sineEigenvectors.shape.rows)).reversedRows()
        let transformedSineEigenvectors = (decomposition <*> sineEigenvectors.padded(top: cosineEigenvectors.shape.rows)).reversedRows()
        
        let eigenvectors = Matrix.init(shape: .square(length: length)) { row, column in
            if length % 2 == 0 && column == length - 1 {
                return transformedCosineEigenvectors[row, column / 2 + 1]
            }
            return column % 2 == 0 ? transformedCosineEigenvectors[row, column / 2] : transformedSineEigenvectors[row, column / 2]
        }
        
        return ComplexMatrix(real: eigenvectors)
    }
    
    public func eigenvalues() -> ComplexMatrix<Scalar> {
        let alpha = order * (.pi / 2)
        let indices = Matrix(shape: .row(length: length)) { row, column in
            return Scalar(column < length - 1 ? column : length - length % 2)
        }
        return ComplexMatrix(real: (indices * -alpha).cos(), imaginary: (indices * -alpha).sin())
    }
    
}

extension DFrFT where Scalar == Double {
    
    internal func hamiltonian(derivativeOrder: Int = 2) -> Matrix<Scalar> {
        let order = derivativeOrder / 2
        var derivative: Matrix<Scalar> = (1...order).reduce(Matrix.zeros(shape: .row(length: length))) { derivative, derivativeOrder in
            var coefficients = finiteDifferenceCoefficients(derivativeOrder: derivativeOrder * 2)
            coefficients = coefficients.padded(right: derivative.shape.columns - coefficients.shape.columns)
            coefficients = coefficients.shifted(columns: -derivativeOrder)
            return derivative + coefficients / (Scalar(derivativeOrder * derivativeOrder) * -coefficients[0, 0] / 2.0)
        }
        derivative[0, 0] = 0.0
        return .circulant(vector: derivative.elements) + .diagonal(vector: derivative.fft1D().real.elements)
    }
    
}

extension DFrFT where Scalar == Double {
    
    internal func finiteDifferenceCoefficients(derivativeOrder: Int) -> Matrix<Scalar> {
        return .init(shape: .row(length: derivativeOrder + 1)) { row, column in
            let index = column
            let sign = index % 2 == 0 ? 1.0 : -1.0
            let coefficient = tgamma(Scalar(derivativeOrder + 1)) / (tgamma(Scalar(index + 1)) * tgamma(Scalar(derivativeOrder - index + 1)))
            return sign * coefficient
        }
    }
    
    internal func decomposition() -> Matrix<Scalar> {
        return .init(shape: .square(length: length)) { row, column in
            let diagonal = column - row
            let antidiagonal = column + row
            switch (diagonal, antidiagonal) {
            case (0, 0):
                return 1.0
            case (0, length):
                return 1.0
            case (0, _) where column + row < length:
                return 1.0 / Scalar.sqrt(2.0)
            case (0, _):
                return -1.0 / Scalar.sqrt(2.0)
            case (_, length):
                return 1.0 / Scalar.sqrt(2.0)
            default:
                return 0.0
            }
        }
    }
    
}
