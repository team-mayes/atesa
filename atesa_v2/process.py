"""
This portion of the program is responsible for handling setup of the appropriate batch script(s) for the next step in a
Thread, passing them to a task manager to submit them, and updating the list of currently running threads accordingly.
"""

def process(thread, running, settings):
    print(thread.coordinates)
