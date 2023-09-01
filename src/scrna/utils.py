import logging

import scrna.cfg as cfg


def get_logger(name, filename):
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(asctime)s | %(message)s")
    ch = logging.StreamHandler()
    fh = logging.FileHandler(cfg.get_root_dir() / "logs" / filename)

    for h in (ch, fh):
        h.setFormatter(formatter)
        logger.addHandler(h)

    return logger
