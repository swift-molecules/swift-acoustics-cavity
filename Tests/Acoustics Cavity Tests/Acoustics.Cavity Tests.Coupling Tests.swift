import Acoustics_Cavity
import Testing

extension `Acoustics.Cavity Tests` {
    @Suite struct Coupling {}
}

extension `Acoustics.Cavity Tests`.Coupling {
    /// The box 0.4 × 0.3 × 0.2 m, closed, with its lid z = 0.2 a piston: one mode moving the lid up by one.
    static let box = `Acoustics.Cavity Tests`.Solid.box(0.4, 0.3, 0.2, cells: (4, 3, 2))

    static var lid: [[Int]] {
        // The top layer's cells' upper faces, corners 4–7.
        box.cells.filter { cell in cell.allSatisfy { box.vertices[$0].z > 0.05 } }.map { Array($0[4..<8]) }
    }

    static func piston(frequency: Double, mass: Double, loss: Double = 0) -> Acoustics.Cavity.Wall {
        let displacement = box.vertices.map { Vector<3, Double>(x: 0, y: 0, z: $0.z == 0.2 ? 1.0 : 0) }
        return Acoustics.Cavity.Wall(
            faces: lid,
            modes: [
                Acoustics.Cavity.Wall.Mode(
                    frequency: Frequency.hertz(frequency),
                    mass: Mass.kilograms(mass),
                    loss: loss,
                    displacement: displacement
                )
            ]
        )
    }

    static let hole = Acoustics.Cavity.Opening(
        area: Area<Space.Physical>.amount(0.006),
        length: Measure<1, Space.Physical>.amount(0.08)
    )

    /// Christensen and Vistisen's two oscillators (J. Acoust. Soc. Am. 68 (1980) 758, Eqs. 2, 5, 7 and 8): a
    /// plate of area A, mass m_p and stiffness k_p and an air piston of area S and mass m_a on a cavity of volume V.
    static func twoOscillators(
        fp0: Double,
        mp: Double,
        area a: Double,
        s: Double,
        ma: Double,
        volume v: Double
    ) -> (Double, Double) {
        let air = Acoustics.Air.room
        let mu = air.density.underlying * air.speed.underlying * air.speed.underlying / v
        let kp = mp * (2 * Double.pi * fp0) * (2 * Double.pi * fp0)
        let p = (kp + mu * a * a) / mp
        let h = mu * s * s / ma
        let ph = mu * mu * a * a * s * s / (mp * ma)
        let root = ((p - h) * (p - h) + 4 * ph).squareRoot()
        return (
            ((p + h - root) / 2).squareRoot() / (2 * Double.pi), ((p + h + root) / 2).squareRoot() / (2 * Double.pi)
        )
    }

    @Test
    func `A piston and a hole on the uniform pressure alone are Christensen and Vistisen's two oscillators`() throws {
        let air = try Acoustics.Cavity.Modes(Self.box, count: 1)
        let body = try Acoustics.Cavity.Coupling(
            air,
            walls: [Self.piston(frequency: 180, mass: 0.06)],
            openings: [Self.hole],
            count: 2
        )
        let ma = 1.204 * 0.006 * 0.08
        let (low, high) = Self.twoOscillators(fp0: 180, mp: 0.06, area: 0.12, s: 0.006, ma: ma, volume: 0.024)
        let found = body.map(\.frequency.underlying)
        #expect(abs(found[0] / low - 1) < 1e-9, "\(found[0]) against \(low)")
        #expect(abs(found[1] / high - 1) < 1e-9, "\(found[1]) against \(high)")
        // The lower mode breathes: the plate and the plug move the same volume out of phase, nearly.
        #expect(body[0].walls[0][0] * body[0].openings[0] < 0)
        #expect(body[1].walls[0][0] * body[1].openings[0] > 0)
    }

    @Test
    func `The squares of the coupled frequencies sum to those of the plate and the Helmholtz resonance`() throws {
        // Christensen and Vistisen's Eq. 9, with the plate's frequency taken with the hole plugged and the
        // Helmholtz frequency with the plate held, each computed here as a coupling of its own.
        let air = try Acoustics.Cavity.Modes(Self.box, count: 1)
        let wall = Self.piston(frequency: 150, mass: 0.05)
        let body = try Acoustics.Cavity.Coupling(air, walls: [wall], openings: [Self.hole], count: 2)
        let plugged = try Acoustics.Cavity.Coupling(air, walls: [wall], count: 1)
        let helmholtz = try Acoustics.Cavity.Coupling(air, walls: [], openings: [Self.hole], count: 1)
        let f = body.map(\.frequency.underlying)
        let fp = plugged[0].frequency.underlying
        let fh = helmholtz[0].frequency.underlying
        #expect(abs((f[0] * f[0] + f[1] * f[1]) / (fp * fp + fh * fh) - 1) < 1e-9)
        #expect(f[0] < Swift.min(fp, fh) && f[1] > Swift.max(fp, fh))
        // A plug alone on the cavity is Helmholtz's resonator: c √(S / V L) / 2π.
        let c = Acoustics.Air.room.speed.underlying
        #expect(abs(fh / (c * (0.006 / (0.024 * 0.08)).squareRoot() / (2 * Double.pi)) - 1) < 1e-9)
    }

    @Test
    func `The admittance has its zero between the two resonances at the Helmholtz frequency`() throws {
        // Christensen and Vistisen's Eq. 6a: the plate's velocity vanishes where the plug alone would ring.
        let air = try Acoustics.Cavity.Modes(Self.box, count: 1)
        let body = try Acoustics.Cavity.Coupling(
            air,
            walls: [Self.piston(frequency: 150, mass: 0.05)],
            openings: [Self.hole],
            count: 2
        )
        let reach = body.reach(at: Self.box.vertices.firstIndex { $0.z == 0.2 }!, along: Vector(x: 0, y: 0, z: 1))
        func compliance(_ f: Double) -> Double {
            var sum = 0.0
            for (mode, r) in zip(body, reach) {
                let w0 = 2 * Double.pi * mode.frequency.underlying
                let w = 2 * Double.pi * f
                sum += r * r / (mode.mass.underlying * (w0 * w0 - w * w))
            }
            return sum
        }
        var (low, high) = (body[0].frequency.underlying * 1.0001, body[1].frequency.underlying * 0.9999)
        #expect(compliance(low) * compliance(high) < 0)
        for _ in 0..<200 {
            let middle = (low + high) / 2
            if compliance(low) * compliance(middle) <= 0 { high = middle } else { low = middle }
        }
        let c = Acoustics.Air.room.speed.underlying
        let fh = c * (0.006 / (0.024 * 0.08)).squareRoot() / (2 * Double.pi)
        #expect(abs(low / fh - 1) < 1e-9, "\(low) against \(fh)")
    }

    @Test
    func `The air's other modes add mass to a wall and lower it`() throws {
        // A piston lid on its own: the uniform pressure stiffens it; the other modes, more of them each time, add
        // the air's mass and pull it down, never up. Only the modes along z feel a uniform lid.
        let wall = Self.piston(frequency: 150, mass: 0.05)
        var previous = Double.infinity
        for count in [1, 4, 8, 16] {
            let air = try Acoustics.Cavity.Modes(Self.box, count: count)
            let body = try Acoustics.Cavity.Coupling(air, walls: [wall], count: 1)
            let f = body[0].frequency.underlying
            #expect(f <= previous * (1 + 1e-12))
            previous = f
        }
        let air = try Acoustics.Cavity.Modes(Self.box, count: 1)
        let stiff = try Acoustics.Cavity.Coupling(air, walls: [wall], count: 1)[0].frequency.underlying
        // Christensen and Vistisen's Eq. 4: f_p² = f_p0² + f_a², f_a² = ρ c² A² / (V m_p) / 4π².
        let rho = Acoustics.Air.room.density.underlying
        let c = Acoustics.Air.room.speed.underlying
        let fa2 = rho * c * c * 0.12 * 0.12 / (0.024 * 0.05) / (4 * Double.pi * Double.pi)
        #expect(abs(stiff * stiff / (150 * 150 + fa2) - 1) < 1e-9)
        // A lid on a shallow box gains little mass, so the drop stays a few per cent.
        #expect(previous < stiff * (1 - 1e-6) && previous / stiff > 0.9, "\(previous) against \(stiff)")
    }

    @Test
    func `A mode loses what its walls lose by their strain energy, and what it radiates`() throws {
        let air = try Acoustics.Cavity.Modes(Self.box, count: 1)
        let body = try Acoustics.Cavity.Coupling(
            air,
            walls: [Self.piston(frequency: 150, mass: 0.05, loss: 0.02)],
            count: 1
        )
        let mode = body[0]
        let rho = Acoustics.Air.room.density.underlying
        let c = Acoustics.Air.room.speed.underlying
        let w = 2 * Double.pi * mode.frequency.underlying
        let kp = 0.05 * (2 * Double.pi * 150) * (2 * Double.pi * 150)
        // One amplitude: the mass is the piston's and the volume its area.
        #expect(abs(mode.mass.underlying - 0.05) < 1e-12)
        #expect(abs(mode.volume.underlying - 0.12) < 1e-12)
        let expected = 0.02 * kp / (w * w * 0.05) + rho * w * 0.12 * 0.12 / (4 * Double.pi * 0.05 * c)
        #expect(abs(mode.loss / expected - 1) < 1e-9)
    }

    @Test
    func `A coupling refuses what it cannot build`() throws {
        let air = try Acoustics.Cavity.Modes(Self.box, count: 3)
        let wall = Self.piston(frequency: 150, mass: 0.05)
        #expect(throws: Acoustics.Cavity.Error.count) {
            try Acoustics.Cavity.Coupling(air, walls: [wall], count: 4)
        }
        let pentagon = Acoustics.Cavity.Wall(faces: [[0, 1, 2, 3, 4]], modes: wall.modes)
        #expect(throws: Acoustics.Cavity.Error.wall(0)) {
            try Acoustics.Cavity.Coupling(air, walls: [pentagon], count: 1)
        }
        let short = Acoustics.Cavity.Wall(
            faces: wall.faces,
            modes: [
                Acoustics.Cavity.Wall.Mode(
                    frequency: Frequency.hertz(150),
                    mass: Mass.kilograms(0.05),
                    loss: 0,
                    displacement: [Vector(x: 0, y: 0, z: 1)]
                )
            ]
        )
        #expect(throws: Acoustics.Cavity.Error.wall(0)) {
            try Acoustics.Cavity.Coupling(air, walls: [short], count: 1)
        }
        // A cavity open somewhere has no uniform pressure for a hole to work against.
        let vented = Acoustics.Cavity(vertices: Self.box.vertices, cells: Self.box.cells, open: [0])
        let open = try Acoustics.Cavity.Modes(vented, count: 2)
        #expect(throws: Acoustics.Cavity.Error.opening) {
            try Acoustics.Cavity.Coupling(open, walls: [], openings: [Self.hole], count: 1)
        }
    }
}
