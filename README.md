# aer302_spring2021_examples_python

## Getting Started
These examples are written in python, so you must install a valid python interpreter alongside the `pip` package manager for python. 

| Linux | MacOS | Windows |
|-------|-----|---------|
|`sudo apt-get install python3 python3-pip `|`brew install python3`|[Download Python here](https://www.python.org/downloads/windows/)

---

A python *virtual environment* is used to run the examples associated with aer302 in an isolated environment. Firstly, install the `virtualenv` environment manager by
```
pip3 install virtualenv
```

Next, create a virtual environment
```
virtualenv aer_env
```
then activate the environment by running
| Linux / MacOS | Windows |
|-|-|
|`source ./aer_env/bin/activate` | `.\aer_env\Scripts\activate` |

---

Lastly, install the dependencies for these examples by running
```
pip3 install -r requirements.txt
```

---

At this step, you should be ready to run the examples! 
| Linux / MacOS | Windows |
|-|-|
| `python3 paperplane.py` | `py paperplane.py` |

---

To run the the performance analysis, you need to specify the plane, altitude, velocity, and the type of analysis (range, endurance, power, climb rate, etc). For example, to compute power required and available
| `python3 performance/examples.py -plane CP1 -H 1000 -V 60 -power`|
Some analysis require an array of altitudes and velocities, see comments in performance/examples.py for detail