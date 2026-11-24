public import Eigen

extension Acoustics.Cavity {
    /// Why a cavity's modes could not be found.
    public enum Error: Swift.Error, Hashable, Sendable {
        /// The cell at this index has neither eight nor four corners, names a vertex that does not exist, or turns
        /// inside out somewhere.
        case cell(Int)
        /// The eigen solver failed.
        case eigen(Eigen.Error)
        /// More modes were asked for than the cavity has.
        case count
    }
}
