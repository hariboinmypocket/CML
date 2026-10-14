#!/bin/sh
# Install this repository's git hooks into the shared hooks directory.
#
# Hooks live in the common git dir, not per worktree, so one install covers the
# main checkout and every worktree under .claude/worktrees/.
#
# The hook is COPIED, not symlinked. Installing from a worktree would otherwise
# leave the hook pointing inside it, and a worktree is temporary -- once it is
# cleaned up git finds a dangling symlink, skips the hook without a word, and
# pushes quietly stop being checked or mirrored. Re-run this script after editing
# scripts/hooks/ to pick the change up.
#
# Uninstall: rm "$(git rev-parse --git-common-dir)/hooks/pre-push"
set -eu

SOURCE_DIR="$(cd "$(dirname "$0")" && pwd)"
HOOK_DIR="$(git rev-parse --git-common-dir)/hooks"

mkdir -p "$HOOK_DIR"
for hook in pre-push; do
    target="$HOOK_DIR/$hook"
    if [ -e "$target" ] && ! grep -q 'installed by scripts/hooks/install.sh' "$target" 2>/dev/null; then
        echo "refusing to replace a hook this script did not install: $target" >&2
        echo "move it aside and re-run." >&2
        exit 1
    fi
    rm -f "$target"
    {
        printf '%s\n' "#!/bin/sh"
        printf '%s\n' "# installed by scripts/hooks/install.sh -- edit scripts/hooks/$hook and re-run."
        tail -n +2 "$SOURCE_DIR/$hook"
    } > "$target"
    chmod +x "$target"
    echo "installed $hook (copied from $SOURCE_DIR/$hook)"
done
