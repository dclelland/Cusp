//
//  DFrFT1D.swift
//  Cusp
//
//  Created by Daniel Clelland on 09/03/2025.
//

import Foundation
import Numerics
import Plinth

extension Matrix where Scalar == Double {
    
    public func dfrft1D(order: Scalar, approximationOrder: Int = 2) -> ComplexMatrix<Scalar> {
        return ComplexMatrix(real: self).dfrft1D(order: order, approximationOrder: approximationOrder)
    }
    
    public func dfrft1D(matrix: ComplexMatrix<Scalar>) -> ComplexMatrix<Scalar> {
        return ComplexMatrix(real: self).dfrft1D(matrix: matrix)
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    public func dfrft1D(order: Scalar, approximationOrder: Int = 2) -> ComplexMatrix<Scalar> {
        let matrix = ComplexMatrix.dfrftMatrix(length: shape.count, order: order, approximationOrder: approximationOrder)
        return dfrft1D(matrix: matrix)
    }
    
    public func dfrft1D(matrix: ComplexMatrix<Scalar>) -> ComplexMatrix<Scalar> {
        return (matrix <*> asColumn()).asRow()
    }
    
}

extension ComplexMatrix where Scalar == Double {
    
    public static func dfrftMatrix(length: Int, order: Scalar, approximationOrder: Int = 2) -> ComplexMatrix<Scalar> {
        let eigenvectors = ComplexMatrix.dfrftEigenvectors(length: length, approximationOrder: approximationOrder)
        let eigenvalues = ComplexMatrix.dfrftEigenvalues(length: length, order: order)
        return eigenvectors <*> .diagonal(vector: eigenvalues.elements) <*> eigenvectors.conjugate().transposed()
    }
    
    public static func dfrftEigenvectors(length: Int, approximationOrder: Int = 2) -> ComplexMatrix<Scalar> {
        let hamiltonian = Matrix.hamiltonian(length: length, approximationOrder: approximationOrder)
        let decomposition = Matrix.decomposition(length: length)
        
        let diagonalized: Matrix = decomposition <*> hamiltonian <*> decomposition.transposed()
        
        let cosines: Matrix = diagonalized[0..<(length / 2 + 1), 0..<(length / 2 + 1)]
        let sines: Matrix = diagonalized[(length / 2 + 1)..<length, (length / 2 + 1)..<length]
        
        let cosineEigenvectors = try! cosines.eigendecomposition(computing: .rightEigenvectors).sorted(.realAscending).rightEigenvectors.real
        let sineEigenvectors = try! sines.eigendecomposition(computing: .rightEigenvectors).sorted(.realAscending).rightEigenvectors.real
        
        let oddCount = Int(ceil(Scalar(length) / 2.0 - 1.0))
        let evenCount = length / 2 + 1
        
        var expandedCosineEigenvectors = Matrix.zeros(shape: .init(rows: oddCount + evenCount, columns: evenCount))
        var expandedSineEigenvectors = Matrix.zeros(shape: .init(rows: oddCount + evenCount, columns: oddCount))
        
        expandedCosineEigenvectors[0..<evenCount, 0..<evenCount] = cosineEigenvectors
        expandedSineEigenvectors[evenCount..<(evenCount + oddCount), 0..<oddCount] = sineEigenvectors
        
        let transformedCosineEigenvectors = (decomposition <*> expandedCosineEigenvectors).reversedRows()
        let transformedSineEigenvectors = (decomposition <*> expandedSineEigenvectors).reversedRows()
        
        let eigenvectors = Matrix.init(shape: .square(length: length)) { row, column in
            if length % 2 == 0 && column == length - 1 {
                return transformedCosineEigenvectors[row, column / 2 + 1]
            }
            return column % 2 == 0 ? transformedCosineEigenvectors[row, column / 2] : transformedSineEigenvectors[row, column / 2]
        }
        
        return ComplexMatrix(real: eigenvectors)
    }
    
    public static func dfrftEigenvalues(length: Int, order: Scalar) -> ComplexMatrix<Scalar> {
        let alpha = order * (.pi / 2)
        let indices = Matrix(shape: .row(length: length)) { row, column in
            return Scalar(column < length - 1 ? column : length - length % 2)
        }
        return ComplexMatrix(real: (indices * -alpha).cos(), imaginary: (indices * -alpha).sin())
    }
    
}

extension Matrix where Scalar == Double {
    
    public static func hamiltonian(length: Int, approximationOrder: Int = 2) -> Matrix<Scalar> {
        let order = approximationOrder / 2
        var sum = (1...order).reduce(Matrix.zeros(shape: .row(length: length))) { sum, k in
            var stencil = Matrix<Scalar>.centralDifference(order: k * 2)
            stencil = stencil.padded(right: sum.shape.columns - stencil.shape.columns)
            stencil = stencil.shifted(columns: -k)
            return sum + stencil / (Scalar(k * k) * -stencil[0, 0] / 2.0)
        }
        sum[0, 0] = 0.0
        return .circulant(vector: sum.elements) + .diagonal(vector: sum.fft1D().real.elements)
    }
    
    public static func centralDifference(order: Int) -> Matrix {
        return .init(shape: .row(length: order + 1)) { row, column in
            let n = order
            let k = column
            return (k % 2 == 0 ? 1.0 : -1.0) * tgamma(Scalar(n + 1)) / (tgamma(Scalar(k + 1)) * tgamma(Scalar(n - k + 1)))
        }
    }
    
    public static func decomposition(length: Int) -> Matrix<Scalar> {
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
