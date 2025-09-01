-- Build an ordered (<ol>) table of contents from document headers.
local function build_ol_toc(blocks)
  -- collect headers in order
  local hs = {}
  for _, b in ipairs(blocks) do
    if b.t == "Header" then
      table.insert(hs, {level=b.level, id=b.identifier, content=b.content})
    end
  end
  if #hs == 0 then return nil end

  -- stack of {level, items} where items is an OrderedList content table
  local stack = { {level=0, items={}} }

  local function push_item(h)
    local link = pandoc.Link(h.content, "#" .. h.id)
    local item = { pandoc.Plain(link) }  -- list item blocks
    table.insert(stack[#stack].items, { item }) -- placeholder list item
  end

  local function start_sublist(level)
    local sub = { }
    -- attach new sublist to the last item of current list
    local cur_items = stack[#stack].items
    local last = cur_items[#cur_items]
    -- add sublist block to that last item
    table.insert(last, pandoc.OrderedList(sub, {start=1, style="Decimal"}))
    table.insert(stack, {level=level, items=sub})
  end

  for _, h in ipairs(hs) do
    local cur = stack[#stack].level
    if h.level > cur then
      -- going deeper
      if #stack[#stack].items == 0 then
        -- first item at this level
        push_item(h)
        start_sublist(h.level)
      else
        start_sublist(h.level)
        push_item(h)
      end
    else
      -- climb up to appropriate level
      while #stack > 1 and stack[#stack].level >= h.level do
        table.remove(stack)
      end
      push_item(h)
    end
  end

  -- top ordered list
  return pandoc.OrderedList(stack[1].items, {start=1, style="Decimal"})
end

-- Replace explicit \tableofcontents with our ordered list TOC.
function TableOfContents(_)  -- triggers when source has \tableofcontents
  local toc = build_ol_toc(PANDOC_DOCUMENT.blocks)
  if toc then
    return pandoc.Div({ toc }, { id = "TOC", class = "toc" })
  end
  return nil
end

