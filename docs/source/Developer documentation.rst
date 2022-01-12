Developer documentation
=======================

Contributing
------------

If you're not familiar with the workflow for contributing code to a GitHub
repository, an excellent place to start is the `AstroPy developer docs
<https://docs.astropy.org/en/stable/development/workflow/development_workflow.html>`_.

It is strongly recommend to use a dedicated development environment if you want to contribute to FULMAR (or any other project).

Set up your development environment
+++++++++++++++++++++++++++++++++++

After getting some familiarity with the workflow, you should fork `the FULMAR
repository <https://github.com/astrojose9/fulmar>`_ and clone it to your
local machine:

.. code-block:: bash

    git clone https://github.com/YOURUSERNAME/fulmar.git
    cd fulmar
    git checkout -b BRANCHNAME

for some name ``BRANCHNAME`` describing your contribution.

Then you should set up an isolated environment for development using a `Conda
environment
<https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_,
`virtualenv <https://virtualenv.pypa.io/>`_, `venv
<https://docs.python.org/3/library/venv.html>`_, or similar.

For `venv <https://docs.python.org/3/library/venv.html>`_ you create your environment executing:

.. code-block:: bash

    python3 -m venv /path/to/new/virtual/environment

For clarity, we will replace ``/path/to/new/virtual/environment`` by ``myenv`` in the following lines.
Remember to activate your environment before using it, and deactivate after:

.. code-block:: bash
    
    # Activate.
    source /path/to/new/virtual/environment/bin/activate
    # Deactivate.
    source /path/to/new/virtual/environment/bin/deactivate

In case you use `Jupyter <https://jupyter.org>`_ in your development workflow, you need to add your development environment to Jupyter.

First, you need to make sure your environment is activated (``source myenv/bin/activate``). 
Next, install `ipykernel <https://github.com/ipython/ipykernel>`_ and add the IPython kernel to Jupyter:

.. code-block:: bash
    
    #install ipykernel
    pip install --user ipykernel

    # Add your virtual environment to Jupyter 
    python -m ipykernel install --user --name=myenv

Congratulations, you can now select your development environment as a kernel in Jupyter.

Testing your contribution
+++++++++++++++++++++++++

**FULMAR** uses `pytest <https://docs.pytest.org/en/stable/>`_ for testing.
Before opening a pull request to the `GitHub repository <https://github.com/astrojose9/Fulmar>`_, you should check if it behaves as expected.
If you're adding new functionality, make sure that you implement at least one test that checks functionality of your pull request. Not only it makes sure your contributions works, but it also prevents further contriutions breaking your cool new feature.

You can run all the `unit tests <https://en.wikipedia.org/wiki/Unit_testing>`_ using the following command:

.. code-block:: bash
    python -m pytest -v tests

Reporting an issue
------------------

If you encounter an issue with Fulmar, we'd be glad if you could `submit an issue on the GitHub repository <https://github.com/astrojose9/fulmar/issues>`_. In case of reporting a bug, please provide a way of reproducing it.