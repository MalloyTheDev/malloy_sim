#include <test_check.hpp>

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <string>
#include <vector>

// Issue #13: the module graph is enforced at LINK time and nowhere else.
//
// Every module declares `target_include_directories(<mod> PUBLIC include)`, so
// the whole `include/` tree is on every consumer's search path. A module can
// therefore include another module's header without linking it, and the build
// succeeds. Nothing catches it, and the rule it breaks is the one the whole
// architecture rests on: domains share only malloy::math and the sim_core
// vocabulary, and no domain knows that any other domain exists (ADR 0006).
//
// This checks the source against the build files, so the two cannot drift.
namespace
{
namespace fs = std::filesystem;

std::string read_file(const fs::path& path)
{
    std::ifstream in(path);
    return std::string((std::istreambuf_iterator<char>(in)),
                       std::istreambuf_iterator<char>());
}

// The text between a `target_link_libraries(` and its matching close paren.
// Both spellings appear in this project: INTERFACE targets put it on one line,
// STATIC ones spread it over several.
std::vector<std::pair<std::string, std::string>> link_blocks(const std::string& text)
{
    std::vector<std::pair<std::string, std::string>> blocks;
    const std::string call = "target_link_libraries(";
    for (std::size_t at = text.find(call); at != std::string::npos;
         at = text.find(call, at + 1))
    {
        std::size_t open = at + call.size();
        std::size_t depth = 1;
        std::size_t scan = open;
        while (scan < text.size() && depth > 0)
        {
            if (text[scan] == '(')
            {
                ++depth;
            }
            else if (text[scan] == ')')
            {
                --depth;
            }
            ++scan;
        }
        const std::string body = text.substr(open, scan - open - 1);

        // The target name is the first token.
        std::size_t end = body.find_first_of(" \t\r\n");
        if (end == std::string::npos)
        {
            continue;
        }
        blocks.emplace_back(body.substr(0, end), body);
    }
    return blocks;
}

// Every `malloy::<name>` in a block, minus the two that carry only flags.
std::set<std::string> aliases_in(const std::string& body)
{
    std::set<std::string> names;
    const std::string prefix = "malloy::";
    for (std::size_t at = body.find(prefix); at != std::string::npos;
         at = body.find(prefix, at + 1))
    {
        std::size_t start = at + prefix.size();
        std::size_t end = start;
        while (end < body.size() &&
               (std::isalnum(static_cast<unsigned char>(body[end])) != 0 ||
                body[end] == '_'))
        {
            ++end;
        }
        const std::string name = body.substr(start, end - start);
        if (name != "project_options" && name != "project_warnings" && !name.empty())
        {
            names.insert(name);
        }
    }
    return names;
}

// Every `#include <malloy/<module>/...>` in one file.
std::set<std::string> malloy_includes(const std::string& text)
{
    std::set<std::string> modules;
    const std::string prefix = "#include <malloy/";
    for (std::size_t at = text.find(prefix); at != std::string::npos;
         at = text.find(prefix, at + 1))
    {
        const std::size_t start = at + prefix.size();
        const std::size_t slash = text.find('/', start);
        if (slash != std::string::npos)
        {
            modules.insert(text.substr(start, slash - start));
        }
    }
    return modules;
}

std::vector<fs::path> sources_of(const fs::path& root, const std::string& module)
{
    std::vector<fs::path> files;
    for (const fs::path& dir :
         {root / "include" / "malloy" / module, root / "src" / module})
    {
        if (!fs::exists(dir))
        {
            continue;
        }
        for (const auto& entry : fs::directory_iterator(dir))
        {
            if (entry.is_regular_file())
            {
                files.push_back(entry.path());
            }
        }
    }
    std::sort(files.begin(), files.end());
    return files;
}
} // namespace

int main()
{
    const fs::path root{MALLOY_REPO_ROOT};
    const std::string cmake = read_file(root / "CMakeLists.txt");
    MALLOY_CHECK_TRUE(!cmake.empty());

    // --- The declared graph, read out of the build files rather than kept in
    //     a table here that could drift from them. ---
    std::map<std::string, std::set<std::string>> declared;
    for (const auto& block : link_blocks(cmake))
    {
        const std::string& target = block.first;
        const std::string prefix = "malloy_";
        if (target.rfind(prefix, 0) != 0 || target.find("_tests") != std::string::npos ||
            target.find("terminal") != std::string::npos)
        {
            continue; // library targets only
        }
        declared[target.substr(prefix.size())] = aliases_in(block.second);
    }
    // A graph this small going empty would make everything below vacuous.
    MALLOY_CHECK_TRUE(declared.size() >= 10);
    MALLOY_CHECK_TRUE(declared.count("charges") == 1);

    // --- Every module's includes must be covered by its declared links.
    //
    //     This is the check the link step cannot make: the whole include tree
    //     is on every consumer's path, so an undeclared include compiles. ---
    for (const auto& entry : declared)
    {
        const std::string& module = entry.first;
        const std::set<std::string>& allowed = entry.second;

        for (const fs::path& file : sources_of(root, module))
        {
            for (const std::string& included : malloy_includes(read_file(file)))
            {
                if (included == module || allowed.count(included) == 1)
                {
                    continue;
                }
                std::cerr << "malloy_" << module << " includes <malloy/" << included
                          << "/...> in " << file.filename().string()
                          << " but does not link malloy::" << included << '\n';
                return 1;
            }
        }
    }

    // --- ADR 0006's central rule, stated directly: no physics domain may
    //     reference another. The check above would allow it if someone added
    //     the link, so this one is about the architecture rather than the
    //     build files agreeing with themselves. ---
    const std::vector<std::string> domains = {"nbody", "particles", "rigid",
                                              "springs", "charges"};
    for (const std::string& module : domains)
    {
        MALLOY_CHECK_TRUE(declared.count(module) == 1);
        for (const std::string& other : domains)
        {
            if (other == module)
            {
                continue;
            }
            if (declared[module].count(other) == 1)
            {
                std::cerr << "malloy_" << module << " links malloy::" << other
                          << ", but no domain may know another exists (ADR 0006)"
                          << '\n';
                return 1;
            }
            for (const fs::path& file : sources_of(root, module))
            {
                if (malloy_includes(read_file(file)).count(other) == 1)
                {
                    std::cerr << "malloy_" << module << " includes <malloy/" << other
                              << "/...> in " << file.filename().string()
                              << ", but no domain may know another exists (ADR 0006)"
                              << '\n';
                    return 1;
                }
            }
        }

        // And each domain really does rest on the shared vocabulary, so the
        // loop above is checking something rather than passing on empty sets.
        MALLOY_CHECK_TRUE(declared[module].count("math") == 1);
        MALLOY_CHECK_TRUE(declared[module].count("sim_core") == 1);
    }

    // --- malloy_sim_core must stay tiny, which CLAUDE.md rule 13 states as
    //     "adding a domain must not grow it". It may depend on math and
    //     nothing else. ---
    MALLOY_CHECK_TRUE(declared.count("sim_core") == 1);
    for (const std::string& dep : declared["sim_core"])
    {
        if (dep != "math")
        {
            std::cerr << "malloy_sim_core depends on malloy::" << dep
                      << "; it may depend on malloy::math and nothing else" << '\n';
            return 1;
        }
    }

    std::cout << "malloy_layering_tests passed\n";
    return 0;
}
