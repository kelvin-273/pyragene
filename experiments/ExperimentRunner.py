from multiprocessing import Process
from flask import Flask
import json


class ExperimentRunner:

    """Multithreaded experiment runner that schedules experiments."""

    def __init__(self, cases, solvers):
        self._cases = list(cases)

        self._n_cases = len(self._cases)

        self._locked = [False] * self._n_cases
        self._finished = [False] * self._n_cases
        self._results = [None] * self._n_cases
        self._complete = False

        self._save_file = None
        self._server = Flask(__name__)

        @self._server.route("/pause")
        def _pause():
            return "wtf"

        @self._server.route("/continue")
        def _continue():
            return "uh-oh"

        @self._server.route("/results")
        def _results():
            return json.dumps(self.current_results())

    def run_server(self):
        self._server.run(debug=True)

    def current_results(self):
        return [
            (self._cases[i], self._results[i])
            for i in range(self._n_cases)
            if self._finished[i]
        ]

    def next_index(self):
        pass

    def pause(self):
        pass

    def play(self):
        pass

    def _run_experiments(self):
        pass



if __name__ == '__main__':
    from time import sleep
    runner = ExperimentRunner(range(100), [
        lambda x: (lambda: (sleep(), 2*x))()[1]
    ])
    runner.run_server()
