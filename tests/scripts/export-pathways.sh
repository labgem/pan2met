#!/usr/bin/env bash
org_id="${1}"
file="${2}"
pathway-tools -lisp -eval "(progn
	(load \"scripts/export-pathways.lisp\")
	(select-organism :org-id '${org_id})
	(write-to-file \"${file}\"
		(format nil \"~{~A~^~%~}\"
			(loop for p in (all-pathways)
				collect (get-frame-name p))))
	(exit))"
