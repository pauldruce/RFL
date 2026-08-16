#!/bin/bash
# Read the incoming tool call JSON payload from stdin
PAYLOAD=$(cat)

# Extract tool details using jq
TOOL_NAME=$(echo "$PAYLOAD" | jq -r '.toolCall.name')

# 1. Check absolute paths for file tools
if [[ "$TOOL_NAME" == "write_to_file" || "$TOOL_NAME" == "replace_file_content" || "$TOOL_NAME" == "delete_file" ]]; then
  TARGET=$(echo "$PAYLOAD" | jq -r '.toolCall.args.TargetFile // .toolCall.args.AbsolutePath')
  if [[ "$TARGET" == /Users/paul/Dev/RFL/research* ]]; then
    echo '{"decision": "deny", "reason": "SECURITY BLOCK: Modifying files in the /research directory is strictly forbidden."}'
    exit 0
  fi
fi

# 2. Check run_command for Cwd and relative paths
if [[ "$TOOL_NAME" == "run_command" ]]; then
  CWD=$(echo "$PAYLOAD" | jq -r '.toolCall.args.Cwd')
  CMD=$(echo "$PAYLOAD" | jq -r '.toolCall.args.CommandLine')
  
  # Block if they are running the command from inside the research folder
  if [[ "$CWD" == /Users/paul/Dev/RFL/research* ]]; then
    echo '{"decision": "deny", "reason": "SECURITY BLOCK: Running commands inside the /research directory is forbidden."}'
    exit 0
  fi
  
  # Block if the command line references the research folder via relative or absolute paths
  # Matches: /Users/paul/Dev/RFL/research, ./research, ../research, research/
  if echo "$CMD" | grep -qE "(/Users/paul/Dev/RFL/research|(^|[[:space:]])(\.\.?/)+research|(^|[[:space:]])research/)"; then
    echo '{"decision": "deny", "reason": "SECURITY BLOCK: The command references the /research directory which is protected."}'
    exit 0
  fi
fi

# Allow everything else
echo '{"decision": "allow"}'
