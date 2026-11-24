public import Complex
public import Galerkin
internal import Numeric

extension Acoustics.Cavity.Coupling {
    /// What each mode sends far away toward `direction`, a unit vector in the cavity's frame, as a volume: the
    /// volume a simple source at the origin would push out to send the same.
    ///
    /// Every part of the walls is a simple source of its own outward motion, and each opening's plug one at its
    /// centre, each phased by where it is: V(û) = Σ ∫ (u · n) e^{−ik û·y} dS + Σ S x e^{−ik û·c}, with
    /// k = ω / c at the mode's own frequency and time going as e^{−iωt}. Moving with acceleration ä per unit of
    /// its amplitude, the mode sends ρ V(û) ä / (4π r) at a distance r far away, less the delay r / c.
    ///
    /// At low frequency V(û) is the mode's `volume`, the simple source the coupling's radiation loss is taken from.
    /// The pressure the air puts back on the walls, which a boundary-element solve adds, is left out: a large flat
    /// part well above its coincidence sends twice as much as it gives here (Rayleigh's baffled board).
    public func radiated(
        toward direction: Vector<3, Double>
    ) -> [Volume<Space.Physical>.Value<Complex.Number<Double>>] {
        let points = air.cavity.vertices
        let speed = air.air.speed.underlying
        let square = Galerkin.Reference.Quadrilateral.quadrature(degree: 5, Double.self)
        let triangle = Galerkin.Reference.Triangle.quadrature(degree: 4, Double.self)
        let squareShapes = square.points.map { Galerkin.Reference.Quadrilateral.shape(at: $0) }
        let triangleShapes = triangle.points.map { Galerkin.Reference.Triangle.shape(at: $0) }
        // Per wall, each face's points: where each lies along `direction`, its normal and area weight, and the
        // shape functions there.
        typealias Sample = (along: Double, area: Vector<3, Double>, values: [Double])
        let samples: [[[Sample]]] = walls.map { wall in
            wall.faces.map { face in
                let shapes = face.count == 4 ? squareShapes : triangleShapes
                let weights = face.count == 4 ? square.weights : triangle.weights
                let corners = face.map { points[$0] }
                return zip(shapes, weights).map { shape, weight in
                    var place = Vector<3, Double>.zero
                    for (value, corner) in zip(shape.values, corners) { place += corner.coordinates.scaled(by: value) }
                    let area = Galerkin.Jacobian.Surface(shape, corners: corners).normal.scaled(by: weight)
                    return (place.dot(direction), area, shape.values)
                }
            }
        }
        return elements.map { mode in
            let k = 2 * Double.pi * mode.frequency.underlying / speed
            var real = 0.0
            var imaginary = 0.0
            for (index, wall) in walls.enumerated() {
                // The mode's motion at the wall's vertices: its wall modes' shapes, weighted by their shares.
                var motion: [Int: Vector<3, Double>] = [:]
                for face in wall.faces {
                    for vertex in face where motion[vertex] == nil {
                        var sum = Vector<3, Double>.zero
                        for (share, wallMode) in zip(mode.walls[index], wall.modes) where share != 0 {
                            sum += wallMode.displacement[vertex].scaled(by: share)
                        }
                        motion[vertex] = sum
                    }
                }
                for (face, points) in zip(wall.faces, samples[index]) {
                    for sample in points {
                        var moved = Vector<3, Double>.zero
                        for (value, vertex) in zip(sample.values, face) {
                            moved += (motion[vertex] ?? .zero).scaled(by: value)
                        }
                        let pushed = moved.dot(sample.area)
                        let phase = -k * sample.along
                        real += pushed * Double.math.cos(phase)
                        imaginary += pushed * Double.math.sin(phase)
                    }
                }
            }
            for (opening, plug) in zip(openings, mode.openings) {
                let pushed = opening.area.underlying * plug
                let phase = -k * opening.centre.coordinates.dot(direction)
                real += pushed * Double.math.cos(phase)
                imaginary += pushed * Double.math.sin(phase)
            }
            return Volume<Space.Physical>.amount(Complex.Number(real, imaginary))
        }
    }
}
