"""Guard against whitespace-typo Dash callback props (the 'value ' bug class).

`Output("thresh_drop", "value ")` (trailing space) is accepted by Dash at import
but the callback output silently never reaches the real `value` prop -- the reset
button never cleared the threshold dropdown, with no error. This is the same
"declared-but-dead" class as the IntOGen column: structurally present, silently
inert. We scan the pages/ source statically (no import -- the pages modules need
the Dash app context) and assert every Output/Input/State prop name is whitespace-
clean and non-empty.
"""
import glob
import os
import re
import unittest

_PAGES = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "pages")
# Output("component-id", "prop")  /  Input(...)  /  State(...)  -- capture the prop (2nd str arg)
_RX = re.compile(r'\b(Output|Input|State)\(\s*["\'][^"\']+["\']\s*,\s*["\']([^"\']*)["\']')


class TestDashCallbackProps(unittest.TestCase):
    def test_no_whitespace_or_empty_callback_props(self):
        offenders = []
        for path in sorted(glob.glob(os.path.join(_PAGES, "*.py"))):
            with open(path, encoding="utf-8") as fh:
                for lineno, line in enumerate(fh, 1):
                    for kind, prop in _RX.findall(line):
                        if prop != prop.strip() or prop == "":
                            offenders.append(f"{os.path.basename(path)}:{lineno} {kind}(...,'{prop}')")
        self.assertEqual(
            offenders, [],
            "Dash callback props with stray whitespace / empty name (silently inert): "
            + "; ".join(offenders),
        )


if __name__ == "__main__":
    unittest.main()
