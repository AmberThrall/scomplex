use thiserror::Error;

#[derive(Error, Debug)]
pub enum OrientationError {
    #[error("cannot swap orientation of unoriented simplex")]
    CannotSwapOrientation,
    #[error("cannot induce orientation from non-coface simplex")]
    NotACoface,
    #[error("cannot induce orientation of unoriented simplex")]
    CannotInduceOrientation,
    #[error("complex is unorientable")]
    Unorientable,
}