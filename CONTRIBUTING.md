# Contributing to [Oimodeler]

Thank you for considering contributing to [Oimodeler]!
We welcome contributions from everyone. Please submit them as pull requests on GitHub.

## How Can I Contribute?

### 🐛 Reporting Bugs

Before creating bug reports, please check existing issues to avoid duplicates.
When you create a bug report, include as many details as possible.

**Great bug reports include:**

- A clear, descriptive title
- Steps to reproduce the behavior
- Expected behavior vs actual behavior
- Environment details (OS, python version, etc.)
- Screenshots (if applicable)

### 💡 Suggesting Features

To see whether a feature is a good idea and what is the best way to implement it,
you can open an issue/discussion on GitHub.

**Great feature requests include:**

- Clear problem statement
- Proposed solution
- Alternative solutions you've considered
- Additional context

### 📝 Improving Documentation

Documentation improvements are always welcome! This includes:

- Fixing typos
- Adding examples
- Clarifying confusing sections
- Translating documentation

The documentation style loosely follows the [numpydocs](https://numpydoc.readthedocs.io/en/latest/format.html).

### 🔧 Submitting Code

We welcome additions to the code-base of [Oimodeler].

### Pull-request Checklist

- [ ] My code follows the project's style guidelines
- [ ] I have performed a self-review of my own code
- [ ] I have commented my code, particularly in hard-to-understand areas
- [ ] I have made corresponding changes to the documentation
- [ ] I have added tests that prove my fix is effective or that my feature works
- [ ] New and existing unit tests pass locally with my changes

## Style Guide

### Commit Messages

A good base for commit messages is [Conventional Commits](https://conventionalcommits.org/):

```
<type>[optional scope](): <description>

[optional body]

[optional footer]
```

#### Types

- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation only
- `style`: Formatting, missing semicolons, etc.
- `refactor`: Code change that neither fixes a bug nor adds a feature
- `test`: Adding missing tests
- `chore`: Maintenance tasks

#### Examples

```
feat(auth): add OAuth2 support
fix(api): handle null response from payment provider
docs(readme): update installation instructions
```
