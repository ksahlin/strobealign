// Based on example in the log crate documentation

use log::{Record, Level, Metadata, SetLoggerError, LevelFilter};

struct SimpleLogger;

impl log::Log for SimpleLogger {
    fn enabled(&self, metadata: &Metadata) -> bool {
        metadata.level() <= Level::Info
    }

    fn log(&self, record: &Record) {
        if self.enabled(record.metadata()) {
            // Print INFO messages without preceding them with "INFO"
            if record.level() == Level::Info {
                eprintln!("{}", record.args());
            } else {
                eprintln!("{} {}", record.level(), record.args())
            }
        }
    }

    fn flush(&self) {}
}

static LOGGER: SimpleLogger = SimpleLogger;

pub fn init() -> Result<(), SetLoggerError> {
    log::set_logger(&LOGGER)
        .map(|()| log::set_max_level(LevelFilter::Info))
}
