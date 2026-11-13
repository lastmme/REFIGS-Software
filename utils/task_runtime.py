"""Process entry point, durable logs and cooperative cancellation."""
from contextlib import redirect_stdout, redirect_stderr
from datetime import datetime
import traceback


class Cancelled(BaseException):
    pass


class TaskFailed(Exception):
    pass


class TaskLog:
    def __init__(self, stream, events, cancel):
        self.stream, self.events, self.cancel = stream, events, cancel

    def checkpoint(self):
        if self.cancel.is_set():
            raise Cancelled()

    def write(self, message):
        self.checkpoint()
        if message:
            line = f'[{datetime.now():%H:%M:%S}] {message.rstrip()}\n'
            self.stream.write(line)
            self.stream.flush()
            self.events.put({'type': 'log', 'text': line})
        return len(message)

    def flush(self):
        self.stream.flush()

    def put(self, message):
        self.write(str(message))


def execute_task(function, params, path, events, cancel):
    status, detail = 'success', ''
    try:
        with open(path, 'a', encoding='utf-8', buffering=1) as stream:
            log = TaskLog(stream, events, cancel)
            try:
                with redirect_stdout(log), redirect_stderr(log):
                    log.checkpoint()
                    function(**params, text_queue=log)
                    log.checkpoint()
            except Cancelled:
                status, detail = 'cancelled', '用户终止任务'
            except TaskFailed as exc:
                status, detail = 'failed', str(exc)
            except BaseException:
                status, detail = 'failed', traceback.format_exc()
            stream.write(f'\n[{status}] {detail}\n')
            stream.flush()
    except BaseException:
        status, detail = 'failed', traceback.format_exc()
    events.put({'type': 'result', 'status': status, 'detail': detail})
