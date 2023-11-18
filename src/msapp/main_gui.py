from msapp.controller.controller import Controller

def run() -> None:
    print("Starting the msa proteoform profiler GUI...")
    controller = Controller()
    controller.start()
