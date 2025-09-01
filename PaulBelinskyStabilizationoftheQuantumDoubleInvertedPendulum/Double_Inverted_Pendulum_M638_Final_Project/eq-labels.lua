function Math(el)
  if el.mathtype == "DisplayMath" then
    local lbl = el.text:match("\\label{%s*([^}]+)%s*}")
    if lbl then
      -- remove the \label{...} from the math content
      el.text = el.text:gsub("\\label{%s*"..lbl.."%s*}", "")
      -- wrap the math in a span with id=label
      return pandoc.Span({ el }, { id = lbl })
    end
  end
end

