// Based on example in the log crate documentation

use log::{Record, Level, Metadata, SetLoggerError};

struct SimpleLogger {
    level: Level,
}

impl SimpleLogger {
    const fn new(level: Level) -> Self { 
        SimpleLogger { level } 
    }
}

impl log::Log for SimpleLogger {
    fn enabled(&self, metadata: &Metadata) -> bool {
        metadata.level() <= self.level
    }

    fn log(&self, record: &Record) {
        if self.enabled(record.metadata()) {
            eprintln!("{}", record.args())
        }
    }
    
    fn flush(&self) {}
}

pub fn init(level: Level) -> Result<(), SetLoggerError> {
    let logger = SimpleLogger::new(level);
    log::set_boxed_logger(Box::new(logger))
        .map(|()| log::set_max_level(level.to_level_filter()))
}
