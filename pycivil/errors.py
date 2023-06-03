class SectionError(Exception):
    def __init__(self, message) -> None:
        super().__init__(message)
        
        
class RebarCoordsError(Exception):
    def __init__(self, message) -> None:
        super().__init__(message)


class OptimizationError(Exception):
    def __init__(self, message) -> None:
        super().__init__(message)