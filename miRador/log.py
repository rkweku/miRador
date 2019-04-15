import logging
import sys

def setupLogger(name):
    """A function to setup a logger for a method

    Args:
        name: The name of the module being logged

    Returns:
        logger instance

    """

    # Initialize our logger
    formatter = logging.Formatter(fmt="%(asctime)s - %(module)s - %(name)s - "\
        "%(message)s")
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    logger.addHandler(handler)

    # Setup the file handler
    fileHandler = logging.FileHandler('miRador.log')
    fileHandler.setFormatter(formatter)
    logger.addHandler(fileHandler)

    return(logger)

def closeLogger(logger):
    """A function to close the logger for a method

    Args:
        logger: The logger instance

    """

    handlerList = logger.handlers[:]

    for handler in handlerList:
        handler.close()
        logger.removeHandler(handler)
