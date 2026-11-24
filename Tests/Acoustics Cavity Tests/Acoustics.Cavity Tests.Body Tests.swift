import Acoustics_Cavity
import Acoustics_Plate
import Mesh
import Testing

extension `Acoustics.Cavity Tests` {
    /// Whole bodies: a spruce top with its hole, on a rigid box (a guitar's) or a rigid bowl (a lute's), with the
    /// air inside, from their geometry alone.
    @Suite struct Body {}
}

extension `Acoustics.Cavity Tests`.Body {
    typealias Plane = Point<2, Double>

    /// A closed outline through `points`, each span a cubic whose tangents are those of the neighbours.
    static func outline(_ points: [Plane]) -> [Bezier<Plane>] {
        let count = points.count
        func at(_ index: Int) -> Plane { points[(index + count) % count] }
        return (0..<count).map { index in
            let (before, start, end, after) = (at(index - 1), at(index), at(index + 1), at(index + 2))
            return Bezier(
                start: start,
                control1: Point(x: start.x + (end.x - before.x) / 6, y: start.y + (end.y - before.y) / 6),
                control2: Point(x: end.x - (after.x - start.x) / 6, y: end.y - (after.y - start.y) / 6),
                end: end
            )
        }
    }

    /// A soundboard over the first vertices of a body's air: its plate, the faces it covers, and its hole.
    struct Top {
        let plate: Acoustics.Plate
        /// The plate's vertex i is the air's vertex `vertices[i]`.
        let vertices: [Int]
        /// The plate's faces as the air's vertices.
        let faces: [[Int]]
        /// The faces cut away for the hole, as the air's vertices.
        let hole: [[Int]]
        /// The area cut away.
        let area: Double
    }

    /// The top `mesh` (the air's first vertices) clamped all round, with the faces whose centres lie within
    /// `radius` of `centre` cut away for the sound hole, whose rim is free.
    static func top(
        _ mesh: Mesh<Plane>,
        hole centre: Plane,
        radius: Double,
        thickness: Double,
        material: Acoustics.Plate.Material,
        materials: [Int: Acoustics.Plate.Material] = [:]
    ) throws -> Top {
        var kept: [[Int]] = []
        var keptMaterials: [Int: Acoustics.Plate.Material] = [:]
        var cut: [[Int]] = []
        var area = 0.0
        for (index, face) in mesh.faces.enumerated() {
            let corners = face.vertices.map { Int(bitPattern: $0) }
            let points = corners.map { mesh.vertices[$0] }
            let cx = points.map(\.x).reduce(0, +) / 4
            let cy = points.map(\.y).reduce(0, +) / 4
            if (cx - centre.x) * (cx - centre.x) + (cy - centre.y) * (cy - centre.y) < radius * radius {
                cut.append(corners)
                var twice = 0.0
                for i in 0..<4 { twice += points[i].x * points[(i + 1) % 4].y - points[(i + 1) % 4].x * points[i].y }
                area += twice / 2
            } else {
                if let other = materials[index] { keptMaterials[kept.count] = other }
                kept.append(corners)
            }
        }
        let used = Array(Set(kept.flatMap { $0 })).sorted()
        var renumber: [Int: Int] = [:]
        for (new, old) in used.enumerated() { renumber[old] = new }
        let plateMesh = try Mesh<Plane>(
            vertices: used.map { mesh.vertices[$0] },
            faces: kept.map { face in try Polygon(vertices: face.map { Index(Ordinal(UInt(renumber[$0]!))) }) }
        )
        let rim = Set(cut.flatMap { $0 }).intersection(Set(used))
        let plate = Acoustics.Plate(
            mesh: plateMesh,
            thickness: Measure<1, Space.Physical>.amount(thickness),
            material: material,
            materials: keptMaterials,
            edge: .clamped,
            exceptions: Dictionary(uniqueKeysWithValues: rim.map { (renumber[$0]!, Acoustics.Plate.Edge.free) })
        )
        return Top(plate: plate, vertices: used, faces: kept, hole: cut, area: area)
    }

    /// The top's modes as a wall of the air, which lies below it: a deflection up is out of the air.
    static func wall(
        _ modes: Acoustics.Plate.Modes,
        top: Top,
        vertices count: Int,
        loss: Double
    ) -> Acoustics.Cavity.Wall
    {
        Acoustics.Cavity.Wall(
            faces: top.faces,
            modes: modes.map { mode in
                var displacement = [Double](repeating: 0, count: count)
                for (index, vertex) in top.vertices.enumerated() { displacement[vertex] = mode.shape[index] }
                return Acoustics.Cavity.Wall.Mode(
                    frequency: mode.frequency,
                    mass: mode.mass,
                    loss: loss,
                    displacement: displacement
                )
            }
        )
    }

    /// The air of a volume mesh, its cells in swift-mesh's corner order, which is swift-galerkin's.
    static func cavity(_ volume: Mesh<Point<3, Double>>.Volume) -> Acoustics.Cavity {
        Acoustics.Cavity(
            vertices: volume.vertices,
            cells: volume.cells.map { cell in cell.corners.map { Int(bitPattern: $0) } }
        )
    }

    /// ∫ φ dS over `faces` by the corner rule, which is exact for the area a bilinear face sweeps when it is flat.
    static func swept(_ mode: Acoustics.Cavity.Wall.Mode, faces: [[Int]], cavity: Acoustics.Cavity) -> Double {
        var sum = 0.0
        for face in faces {
            let p = face.map { cavity.vertices[$0] }
            // A flat quadrilateral's area is half the cross product of its diagonals; its bilinear mean is the
            // corners' mean only for a parallelogram, which the caller's faces are close to.
            let (ax, ay) = (p[2].x - p[0].x, p[2].y - p[0].y)
            let (bx, by) = (p[3].x - p[1].x, p[3].y - p[1].y)
            let area = abs(ax * by - ay * bx) / 2
            sum += area * face.map { mode.displacement[$0] }.reduce(0, +) / Double(face.count)
        }
        return sum
    }

    /// The top vertex nearest `point`.
    static func nearest(_ point: Plane, in mesh: Mesh<Plane>) -> Int {
        mesh.vertices.indices.min { a, b in
            let (p, q) = (mesh.vertices[a], mesh.vertices[b])
            return (p.x - point.x) * (p.x - point.x) + (p.y - point.y) * (p.y - point.y)
                < (q.x - point.x) * (q.x - point.x) + (q.y - point.y) * (q.y - point.y)
        }!
    }

    /// |Y| = |Σ iω r² / (m (ω₀² (1 + iη) − ω²))| at `frequency`.
    static func admittance(_ body: Acoustics.Cavity.Coupling, reach: [Double], at frequency: Double) -> Double {
        let w = 2 * Double.pi * frequency
        var (re, im) = (0.0, 0.0)
        for (mode, r) in zip(body, reach) {
            let w0 = 2 * Double.pi * mode.frequency.underlying
            let a = w0 * w0 - w * w
            let b = w0 * w0 * mode.loss
            let scale = r * r / mode.mass.underlying / (a * a + b * b)
            // iω (a − i b) = ω b + i ω a
            re += scale * w * b
            im += scale * w * a
        }
        return (re * re + im * im).squareRoot()
    }
}

extension `Acoustics.Cavity Tests`.Body {
    /// `material` with spruce bars of `width`, and of `height` at each point, glued on across the grain at `bars` along x, each smeared
    /// over the faces it crosses in proportion to how much of each face's width it covers: the bar's bending stiffness
    /// E_L I about the plate's midplane added to the face's cross-grain bending, Q₁₂ kept, its shear and mass added.
    static func barred(
        _ mesh: Mesh<Plane>,
        material: Acoustics.Plate.Material,
        thickness t: Double,
        bars: [Double],
        width b: Double,
        height: (Plane) -> Double
    ) -> [Int: Acoustics.Plate.Material] {
        let wood = Acoustics.Plate.Material.Wood.sitkaSpruce
        let el = wood.young.longitudinal.underlying
        let ex = material.young.x.underlying
        let ey = material.young.y.underlying
        let nxy = material.poisson.underlying
        let d = 1 - nxy * nxy * ey / ex
        let (q12, q22) = (nxy * ey / d, ey / d)
        var result: [Int: Acoustics.Plate.Material] = [:]
        for (index, face) in mesh.faces.enumerated() {
            let xs = face.vertices.map { mesh.vertices[Int(bitPattern: $0)].x }
            let (low, high) = (xs.min()!, xs.max()!)
            var cover = 0.0
            for x in bars {
                cover += Swift.max(0, Swift.min(high, x + b / 2) - Swift.max(low, x - b / 2)) / (high - low)
            }
            guard cover > 0 else { continue }
            let ys = face.vertices.map { mesh.vertices[Int(bitPattern: $0)].y }
            let h = height(Point(x: (low + high) / 2, y: ys.reduce(0, +) / 4))
            let offset = (t + h) / 2
            let inertia = b * h * h * h / 12 + b * h * offset * offset
            let q22b = q22 + cover * 12 * el * inertia / (b * t * t * t)
            let nu = q12 / q22b
            let eyb = q22b / (1 + q22b * nu * nu / ex)
            var barred = material
            barred.young.y = Modulus.Young.amount(eyb)
            barred.poisson = Poisson.amount(nu)
            barred.shear.yz = Modulus.Shear.amount(
                material.shear.yz.underlying + cover * wood.shear.longitudinalRadial.underlying * h / t
            )
            barred.density = Density.amount(material.density.underlying + cover * wood.density.underlying * h / t)
            result[index] = barred
        }
        return result
    }
}

extension `Acoustics.Cavity Tests`.Body {
    /// A classical guitar's outline, its length along x as the top's grain runs: 0.49 m long, 0.37 m across the
    /// lower bout, 0.24 m at the waist and 0.28 m across the upper bout.
    static let guitar = outline([
        Point(x: -0.27, y: 0), Point(x: -0.2, y: -0.16), Point(x: -0.09, y: -0.185), Point(x: 0.03, y: -0.12),
        Point(x: 0.12, y: -0.14), Point(x: 0.2, y: -0.11), Point(x: 0.22, y: 0), Point(x: 0.2, y: 0.11),
        Point(x: 0.12, y: 0.14), Point(x: 0.03, y: 0.12), Point(x: -0.09, y: 0.185), Point(x: -0.2, y: 0.16),
    ])

    /// The guitar: its outline 0.095 m deep, a 2.5 mm quartersawn Sitka top with no braces, clamped all round, and
    /// an 88 mm sound hole, Christensen and Vistisen's (1980) diameter, centred 0.09 m towards the neck.
    struct Guitar {
        let mesh: Mesh<Plane>
        let cavity: Acoustics.Cavity
        let volume: Double
        let air: Acoustics.Cavity.Modes
        let wall: Acoustics.Cavity.Wall
        let hole: Acoustics.Cavity.Opening

        init() throws {
            let mesh = try Mesh(inside: `Acoustics.Cavity Tests`.Body.guitar, around: Point(x: -0.02, y: 0), cells: 10)
            let volume = try Mesh<Point<3, Double>>.Volume(below: mesh, layers: 4) { _ in 0.095 }
            let cavity = `Acoustics.Cavity Tests`.Body.cavity(volume)
            let top = try `Acoustics.Cavity Tests`.Body.top(
                mesh,
                hole: Point(x: 0.09, y: 0),
                radius: 0.044,
                thickness: 0.0025,
                material: Acoustics.Plate.Material(.sitkaSpruce, cut: .quarter)
            )
            self.mesh = mesh
            self.cavity = cavity
            self.volume = volume.volumes().reduce(0, +)
            self.air = try Acoustics.Cavity.Modes(cavity, count: 12)
            self.wall = `Acoustics.Cavity Tests`.Body.wall(
                try Acoustics.Plate.Modes(top.plate, count: 10),
                top: top,
                vertices: cavity.vertices.count,
                loss: 0.01
            )
            self.hole = .circle(radius: Space.Physical.metres(0.044), thickness: Space.Physical.metres(0.0025))
        }
    }

    @Test
    func `The guitar's air rings as the air of measured guitars does`() throws {
        let guitar = try Guitar()
        // About 13 litres, Christensen and Vistisen's guitar's.
        #expect(abs(guitar.volume - 0.013) < 0.0005, "\(guitar.volume)")
        // Rossing (ed.), The Science of String Instruments (2010), Table 3.1, p. 25: the air of a Kohno 30 and a
        // Conrad classical guitar with top, back and ribs held: A1 396 and 391, A2 560 and 558, A3 674 and 711,
        // A4 780 and 772 Hz. The box is not theirs, so within 10% of their means.
        let measured = [393.5, 559, 692.5, 776]
        let found = guitar.air.dropFirst().prefix(4).map(\.frequency.underlying)
        for (f, m) in zip(found, measured) { #expect(abs(f / m - 1) < 0.1, "\(f) against \(m)") }
        // Helmholtz's resonance with the top held: Christensen and Vistisen measured 126.5 Hz on theirs and
        // 122–135 Hz on others (p. 760).
        let helmholtz = try Acoustics.Cavity.Coupling(guitar.air, walls: [], openings: [guitar.hole], count: 1)
        let fh = helmholtz[0].frequency.underlying
        #expect(fh > 122 && fh < 135, "\(fh)")
    }

    @Test
    func `The guitar's lowest mode breathes below both the top's and the air's own`() throws {
        // Christensen and Vistisen (1980), p. 762: f₋ lies below f_h and f_p, and there the top and the air in the
        // hole move against each other, so the volume the top pushes out is largely breathed back in.
        let guitar = try Guitar()
        let body = try Acoustics.Cavity.Coupling(guitar.air, walls: [guitar.wall], openings: [guitar.hole], count: 6)
        let plugged = try Acoustics.Cavity.Coupling(guitar.air, walls: [guitar.wall], count: 1)
        let helmholtz = try Acoustics.Cavity.Coupling(guitar.air, walls: [], openings: [guitar.hole], count: 1)
        let low = body[0]
        #expect(low.frequency.underlying < plugged[0].frequency.underlying)
        #expect(low.frequency.underlying < helmholtz[0].frequency.underlying)
        let topVolume = zip(low.walls[0], guitar.wall.modes).map { share, mode in
            share * `Acoustics.Cavity Tests`.Body.swept(mode, faces: guitar.wall.faces, cavity: guitar.cavity)
        }.reduce(0, +)
        let holeVolume = low.openings[0] * guitar.hole.area.underlying
        #expect(topVolume * holeVolume < 0)
        #expect(abs(low.volume.underlying) < abs(topVolume))
    }
}

extension `Acoustics.Cavity Tests`.Body {
    static let lute = outline([
        Point(x: -0.18, y: 0), Point(x: -0.15, y: -0.11), Point(x: -0.06, y: -0.1575), Point(x: 0.05, y: -0.145),
        Point(x: 0.15, y: -0.095), Point(x: 0.225, y: -0.04), Point(x: 0.252, y: 0), Point(x: 0.225, y: 0.04),
        Point(x: 0.15, y: 0.095), Point(x: 0.05, y: 0.145), Point(x: -0.06, y: 0.1575), Point(x: -0.15, y: 0.11),
    ])

    /// A lute of Gamut's Fantasia model's size, 43.2 cm long, 31.5 cm wide and 15.3 cm deep, its bowl swift-mesh's
    /// under the outline. Firth's (1977) four bars at 10, 20, 30 and 40 cm from the tail of a 45 cm body and his
    /// bridge at 11.5 cm (Rossing 2010, Fig. 3.27b, p. 42), scaled to this length. The bars are Schreiner's on the
    /// Jauch lute, 4.2 mm wide and 15.8 mm high in the middle falling to 6 mm at the ribs, on a 2.0 mm quartersawn
    /// Sitka top; the rose is taken as an open hole of 80 mm between the upper two bars. The bowl is rigid.
    struct Lute {
        static let length = 0.432
        static func along(_ centimetres: Double) -> Double { -0.18 + length * centimetres / 45 }

        let mesh: Mesh<Plane>
        let volume: Double
        let air: Acoustics.Cavity.Modes
        let wall: Acoustics.Cavity.Wall
        let rose: Acoustics.Cavity.Opening

        init() throws {
            let centre = Point<2, Double>(x: -0.02, y: 0)
            let outline = `Acoustics.Cavity Tests`.Body.lute
            let mesh = try Mesh(inside: outline, around: centre, cells: 8)
            let volume = try Mesh<Point<3, Double>>.Volume(under: outline, around: centre, depth: 0.153, cells: 8)
            let cavity = `Acoustics.Cavity Tests`.Body.cavity(volume)
            let spruce = Acoustics.Plate.Material(.sitkaSpruce, cut: .quarter)
            let bars = `Acoustics.Cavity Tests`.Body.barred(
                mesh,
                material: spruce,
                thickness: 0.002,
                bars: [10.0, 20, 30, 40].map(Self.along),
                width: 0.0042
            ) { p in
                let half = 0.1575 * Swift.max(0, 1 - (p.x + 0.02) * (p.x + 0.02) / (0.24 * 0.24)).squareRoot()
                let share = Swift.min(1, abs(p.y) / Swift.max(half, 0.01))
                return 0.006 + 0.0098 * (1 - share * share)
            }
            let top = try `Acoustics.Cavity Tests`.Body.top(
                mesh,
                hole: Point(x: Self.along(35), y: 0),
                radius: 0.04,
                thickness: 0.002,
                material: spruce,
                materials: bars
            )
            self.mesh = mesh
            self.volume = volume.volumes().reduce(0, +)
            self.air = try Acoustics.Cavity.Modes(cavity, count: 12)
            self.wall = `Acoustics.Cavity Tests`.Body.wall(
                try Acoustics.Plate.Modes(top.plate, count: 12),
                top: top,
                vertices: cavity.vertices.count,
                loss: 0.01
            )
            self.rose = .circle(radius: Space.Physical.metres(0.04), thickness: Space.Physical.metres(0.002))
        }
    }

    @Test
    func `The lute's lowest resonance at the bridge is its air's, near Firth's`() throws {
        let lute = try Lute()
        // Two thirds of the depth times the outline's area, about 10 litres.
        #expect(abs(lute.volume - 0.0101) < 0.0005, "\(lute.volume)")
        let body = try Acoustics.Cavity.Coupling(lute.air, walls: [lute.wall], openings: [lute.rose], count: 8)
        // The input admittance at the treble end of the bridge, whose lowest peak Firth (1977) measured at 132 Hz
        // and named the air's (Rossing 2010, §3.10.1, p. 41). The lute is not Firth's, so within 15%.
        let bridge = `Acoustics.Cavity Tests`.Body.nearest(Point(x: Lute.along(11.5), y: -0.05), in: lute.mesh)
        let reach = body.reach(at: bridge)
        let frequencies = stride(from: 60.0, through: 700, by: 0.5).map { $0 }
        let magnitude = frequencies.map { `Acoustics.Cavity Tests`.Body.admittance(body, reach: reach, at: $0) }
        let peak = (1..<(frequencies.count - 1)).first {
            magnitude[$0] > magnitude[$0 - 1] && magnitude[$0] > magnitude[$0 + 1]
        }!
        #expect(abs(frequencies[peak] / 132 - 1) < 0.15, "\(frequencies[peak])")
        // It is the air's: the rose breathes most there, and the top little.
        let low = body[0]
        #expect(abs(low.frequency.underlying - frequencies[peak]) < 1)
        #expect(abs(low.openings[0]) == 1)
        // Above it the top's own modes ring, the first one well above the air's.
        #expect(body[1].frequency.underlying > 2 * low.frequency.underlying)
    }
}
