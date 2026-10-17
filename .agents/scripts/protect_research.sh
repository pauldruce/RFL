#!/bin/bash
# Read the incoming tool call JSON payload from stdin
PAYLOAD=$(cat)

# Extract tool details using jq
TOOL_NAME=$(echo "$PAYLOAD" | jq -r '.toolCall.name')

WORKSPACE_ROOT=$(cd "$(dirname "$0")/../.." && pwd)

# 1. Check absolute paths for file tools
if [[ "$TOOL_NAME" == "write_to_file" || "$TOOL_NAME" == "replace_file_content" || "$TOOL_NAME" == "delete_file" ]]; then
  TARGET=$(echo "$PAYLOAD" | jq -r '.toolCall.args.TargetFile // .toolCall.args.AbsolutePath')
  if [[ "$TARGET" == "$WORKSPACE_ROOT/research"* ]]; then
    echo '{"decision": "deny", "reason": "SECURITY BLOCK: Modifying files in the /research directory is strictly forbidden."}'
    exit 0
  fi
fi

# 2. Check run_command for Cwd and relative paths
if [[ "$TOOL_NAME" == "run_command" ]]; then
  CWD=$(echo "$PAYLOAD" | jq -r '.toolCall.args.Cwd')
  CMD=$(echo "$PAYLOAD" | jq -r '.toolCall.args.CommandLine')

  # Block if they are running the command from inside the research folder
  if [[ "$CWD" == "$WORKSPACE_ROOT/research"* ]]; then
    echo '{"decision": "deny", "reason": "SECURITY BLOCK: Running commands inside the /research directory is forbidden."}'
    exit 0
  fi

  # Block if the command line references the research folder via relative or absolute paths
  # Matches: <absolute_workspace_root>/research, ./research, ../research, research/
  if echo "$CMD" | grep -qE "($WORKSPACE_ROOT/research|(^|[[:space:]])(\.\.?/)+research|(^|[[:space:]])research/)"; then
    echo '{"decision": "deny", "reason": "SECURITY BLOCK: The command references the /research directory which is protected."}'
    exit 0
  fi

  # Auto-allow read-only and safe verification mise commands
  if echo "$CMD" | grep -qE "^mise[[:space:]]+(tasks|ls|doctor|version|run[[:space:]]+(affected|plan|ci(:[a-z0-9_:-]+)?|lint(:[a-z0-9_:-]+)?|provision(:[a-z0-9_:-]+)?|traceability:(matrix|query)|test(:[a-z0-9_:-]+)?|pre[:\-]commit|pre[:\-]push|changelog:preview|release:verify))([[:space:]]|$)"; then
    echo '{"decision": "allow", "permissionOverrides": ["command(mise)"], "overwrite": {"BypassSandbox": true}}'
    exit 0
  fi
fi

# Allow everything else
echo '{"decision": "allow"}'

