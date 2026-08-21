
# Notes on Running And Writing Tests

The testing suite is executed using [pytest] with all the tests
inside the `tests/` directory.

## Dependencies

To install the optional `test` dependencies either execute `pip`

```bash
pip install -e .[test]
```

or make use of the [`pixi`](https://pixi.prefix.dev/latest/) framework

```bash
pixi run test-slow
```

## Slow Tests

There are a few slow tests. To exclude the slow tests [pytest] can be run as follows

```bash
pytest -m 'not slow'
```

or

```bash
pixi run test
```

Long tests include running all scripts and notebooks contained
in the `examples/` directory

## Parallelisation

The tests can be run in parallel (e.g. with 10 CPUs) using

```bash
pytest  -m 'not slow' --dist loadgroup  -n 10
```

or

```bash
pixi run test-parallel
```
