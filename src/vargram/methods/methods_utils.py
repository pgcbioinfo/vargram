"""Module to perform auxiliary tasks for VARGRAM class methods."""

class CustomError(Exception):
    """Base class for custom exceptions."""
    pass

class UsageError(CustomError):
    """Raised when a specific condition is not met."""

    def __init__(self, message="Check ordering of method calls."):
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"{self.__class__.__name__}: {self.message}"
    
def usage_checker(bar):
    """Checks whether VARGRAM usage is correct in terms of methods called.
    
    Args:
        bar (bool): True if method bar() has been called.
    Returns:
        bool: True if usage is correct and False otherwise.
    """

    # Return error if bar is not present
    if bar == False:
        return False