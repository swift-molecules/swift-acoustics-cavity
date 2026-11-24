import Acoustics_Cavity
import Testing

extension `Acoustics.Cavity Tests` {
    @Suite struct Helmholtz {}
}

extension `Acoustics.Cavity Tests`.Helmholtz {
    typealias Base = `Acoustics.Cavity Tests`

    @Test
    func `A body with a neck rings at the Helmholtz resonance`() throws {
        // Ingard (1953), On the theory and design of acoustic resonators, J. Acoust. Soc. Am. 25, 1037: f =
        // (c / 2π) √(S / V L′), with the neck's length corrected by 0.85 a at its flanged outer end and by
        // 0.85 a (1 − 1.25 a / R) where it opens into a cylindrical body of radius R. Here f = 91.1 Hz.
        let body = (radius: 0.1, depth: 0.1)
        let neck = (radius: 0.01, length: 0.02)
        let a = neck.radius
        let length = neck.length + 0.85 * a + 0.85 * a * (1 - 1.25 * a / body.radius)
        let volume = Double.pi * body.radius * body.radius * body.depth
        let area = Double.pi * a * a
        let f = Base.c / (2 * Double.pi) * (area / (volume * length)).squareRoot()

        var found: [Double] = []
        for fineness in [1, 2] {
            let cavity = Base.resonator(
                body: body,
                neck: neck,
                outside: (radius: 0.3, height: 0.3),
                fineness: fineness
            )
            let modes = try Acoustics.Cavity.Modes(cavity, count: 2, air: Base.air)
            found.append(modes[0].frequency.underlying)
            // The next mode is the outside air's own, far above.
            #expect(modes[1].frequency.underlying > 5 * f)
        }
        // The rim of the neck is singular, so the finer mesh is within 2%; at second order, the extrapolation to
        // no cell size is within 1%.
        #expect(Base.error(found[1], f) < 0.02, "\(found[1]) Hz, Ingard \(f)")
        let extrapolated = (4 * found[1] - found[0]) / 3
        #expect(Base.error(extrapolated, f) < 0.01, "\(extrapolated) Hz, Ingard \(f)")
    }
}
