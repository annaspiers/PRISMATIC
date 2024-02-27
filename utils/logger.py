import datetime as dt
import logging
import os

from configparser import ConfigParser
# from stat import FILE_ATTRIBUTE_HIDDEN


class Logger(logging.Logger):
    def __init__(self, cfg_path):
        config = ConfigParser()
        with open(os.path.abspath(cfg_path)) as cfg:
            log_cfg = 'LOG'
            config.read_file(cfg)
            if config.has_section(log_cfg):
                log_level = config.get(log_cfg, 'log_level')
                self.log_dir = config.get(log_cfg, 'log_dir')
            try:
                log_level = eval(log_level)
            except Exception:
                err_msg = ('Unable to evaluate config file log level: '
                           f'{log_level}. Using default log '
                           'level INFO.')
                print(err_msg)
                log_level = logging.INFO

            self.log_file_timestamp = dt.datetime.now()
            log_ts_fmt = '%Y%m%d%H%M%S'
            timestamp_str = self.log_file_timestamp.strftime(log_ts_fmt)
            self.log_file_name = \
                f'prismatic_inititalizing_report_{timestamp_str}.log'
            if not os.path.exists(self.log_dir):
                os.makedirs(self.log_dir)
            self.log_path = os.path.join(self.log_dir, self.log_file_name)
            default_format = ('%(asctime)s [%(levelname)s] '
                              '%(name)s - %(message)s')
            default_formatter = logging.Formatter(default_format)

            # initialize root level logger
            logger = logging.getLogger('')
            logger.setLevel(log_level)
            file_handler = logging.FileHandler(self.log_path)
            file_handler.setLevel(log_level)
            file_handler.setFormatter(default_formatter)

            console = logging.StreamHandler()
            console.setLevel(log_level)
            console.setFormatter(default_formatter)

            logger.addHandler(file_handler)
            logger.addHandler(console)

    def getLogger(self, name):
        '''Return the logger from logging library'''
        self._log = logging.getLogger(str(name))
        return self

    def getOriginalLogger(self):
        '''Return original logger'''
        return self._log

    def getName(self):
        '''Return the logger name'''
        return self.getOriginalLogger().name

    def get_log_dir(self):
        return self.log_dir
