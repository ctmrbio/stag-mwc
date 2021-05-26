# This file is part of StaG-mwc
 
class UserMessages():
    """
    Deferred user messages for delayed printout
    """

    def __init__(self):
        self.messages = {
            "info": set(),
            "warning": set(),
        }

    def print_messages(self):
        for level, messages in self.messages.items():
            for message in messages:
                print(level.upper()+":", message)

    def print_info(self):
        for message in self.info:
            print(level.upper()+":", message)

    def print_warnings(self):
        for message in self.warn:
            print(level.upper()+":", message)

    def info(self, message):
        self.messages["info"].add(message)

    def warn(self, message):
        self.messages["warning"].add(message)

        
