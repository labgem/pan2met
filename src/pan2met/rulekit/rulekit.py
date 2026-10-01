"""
A rule kit inspired by https://docs.rs/rule_kit/
"""

import queue
from abc import ABC


class Rule(ABC):
    """
    An abstract class to implement a rule
    """

    def __init__(self):
        pass

    @property
    def priority(self) -> int:
        raise NotImplementedError()

    def evaluate(self, context) -> bool:
        raise NotImplementedError()

    def apply(self, context):
        raise NotImplementedError()


class RuleEngine:
    def __init__(self, rules: list[Rule]):
        self.rules = rules

    def evaluate_all(self, context):
        for rule in self.rules:
            if rule.evaluate(context):
                rule.apply(context)

    def evaluate_first(self, context) -> bool:
        priority_queue = queue.PriorityQueue()
        for rule in self.rules:
            priority_queue.put((rule.priority, rule))

        while not priority_queue.empty():
            _priority, rule = priority_queue.get()
            if rule.evaluate(context):
                rule.apply(context)
                return True
        return False
